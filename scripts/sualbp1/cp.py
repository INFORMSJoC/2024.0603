import argparse
import math
import time
import sys

import docplex.cp.model as cp

import utils_sualbsp1


def solve(
    number_of_tasks,
    cycle_time,
    task_times,
    m_lb, m_ub, all_predecessors, all_followers, earliest, latest,
    forward,
    backward,
    backward_array,
    FiF, PiF, FiB, PiB,
    aijF, bijF, aijB, bijB,
    dij, AiNew,
    distance,
    time_limit=None,
    threads=1,
    verbose=False,
    history=None,
):
    time_start = time.time()

    # 0,...,N-1
    tasks = list(range(0, number_of_tasks))
    stations = list(range(1, m_ub+1))
    setup_actual = []
    # 0,...,N+1
    for i in range(number_of_tasks+2):
        setup_actual.append([])
        # 0,...,N+1
        for j in range(number_of_tasks+2):
            setup_actual[i].append(forward[i][j])
            

    # initialize a CP model
    mdl = cp.CpoModel()

    # use binary variables to represent usage of stations
    usage = {}
    for s in range(m_lb+1):
        usage[s] = 1
    for s in range(m_lb+1, m_ub+1):
        usage[s] = cp.integer_var(0, 1)
    for s1 in range(1, len(stations)):
        for s2 in range(s1+1, len(stations)+1):
            mdl.add(usage[s1] >= usage[s2])


    # use alternative interval variables
    t_to_s = {}
    type_of_s = {}
    setup_of_s = {}
    tasks_of_s = {}
    for s in stations:
        t_to_s[s] = {}
        type_of_s[s] = []
        setup_of_s[s] = []
        tasks_of_s[s] = []
        for i in tasks:
            if earliest[i+1] <= s <= latest[i+1]:
                t_to_s[s][i] = cp.interval_var(name='t{}_to_s{}'.format(i, s), 
                                               optional=True, 
                                               size=task_times[i+1])
                type_of_s[s].append(i)
                tasks_of_s[s].append(i)
        # two dummy interval variables
        # dummy start
        t_to_s[s][number_of_tasks] = \
            cp.interval_var(name='t{}_to_s{}'.format(number_of_tasks, s), size=0)
        type_of_s[s].append(number_of_tasks)
        tasks_of_s[s].append(number_of_tasks)
        # dummy end
        t_to_s[s][number_of_tasks+1] = \
            cp.interval_var(name='t{}_to_s{}'.format(number_of_tasks+1, s), size=0)
        type_of_s[s].append(number_of_tasks+1)
        tasks_of_s[s].append(number_of_tasks+1)
        
    for s in stations:
        for i in range(len(tasks_of_s[s])):
            setup_of_s[s].append([])
            for j in range(len(tasks_of_s[s])):
                setup_of_s[s][-1].append(forward[tasks_of_s[s][i]][tasks_of_s[s][j]])
            

    # Build actual tasks as an alternative between machines
    task_actual = [cp.interval_var(name='actual_t{}'.format(i), 
                                   size=task_times[i+1]) for i in tasks]

    # alternative constraint to synchronize interval variables
    for i in tasks:
        mdl.add(cp.alternative(task_actual[i], [t_to_s[s][i] for s in stations 
                                                if earliest[i+1] <= s <= latest[i+1]], 
                               cardinality = 1))

    # Constrain tasks to no overlap on each station
    sequence = {}
    for s in stations:
        sequence[s] = cp.sequence_var(list(t_to_s[s].values()), name='sq_at_s{}'.format(s), types=type_of_s[s])
        mdl.add(cp.no_overlap(sequence[s], setup_actual, is_direct=True))
        mdl.add(cp.first(sequence[s], t_to_s[s][number_of_tasks]))
        mdl.add(cp.last(sequence[s], t_to_s[s][number_of_tasks+1]))
    
    # link usage with interval variables    
    for s in stations:
        for i in tasks:
            if earliest[i+1] <= s <= latest[i+1]:
                mdl.add(cp.start_of(t_to_s[s][i]) >= 0)
                mdl.add(cp.presence_of(t_to_s[s][i]) <= usage[s])
        mdl.add(cp.start_of(t_to_s[s][number_of_tasks]) >= 0)
        mdl.add(cp.start_of(t_to_s[s][number_of_tasks+1]) >= 0)
        mdl.add(cp.end_of(t_to_s[s][number_of_tasks+1]) + 
            cp.element(backward_array, 
                cp.type_of_prev(sequence[s], t_to_s[s][number_of_tasks+1]) * (number_of_tasks+2) + 
                cp.type_of_next(sequence[s], t_to_s[s][number_of_tasks]) )
            <= 
            cycle_time * usage[s]
        )
   
    # precedence constraints: i precedes j
    for i in tasks:
        for j in all_followers[i+1]:
            for s in stations:
                if earliest[i+1] <= s <= latest[i+1] and earliest[j] <= s <= latest[j]:
                    mdl.add(cp.end_before_start(t_to_s[s][i], t_to_s[s][j-1]))
            mdl.add(sum(s * cp.presence_of(t_to_s[s][i]) for s in stations 
                        if earliest[i+1] <= s <= latest[i+1]) 
                    + distance[i+1, j]
                    <= 
                    sum(s * cp.presence_of(t_to_s[s][j-1]) for s in stations 
                        if earliest[j] <= s <= latest[j]))

    # Minimize the makespan
    mdl.add(cp.minimize(sum([usage[s] for s in stations])))
    
    # Solve model
    is_infeasible = False
    is_optimal = False
    cost = -1
    best_bound = -1
    assignment = []
    res = []

    print('Solving model...')
    solver = cp.CpoSolver(
        mdl, TimeLimit=time_limit, Workers=1
    )
    is_new_solution = True
    while is_new_solution and time.time() - time_start <= time_limit:
        result = solver.search_next()
        is_new_solution = result.is_new_solution()

        if is_new_solution:
            res.append(
                (time.time() - time_start, 
                result.get_objective_value(),
                result.get_objective_bound())
            )

    time_used = time.time() - time_start
    time_out = time.time() - time_start

    if result.is_solution():
        cost = result.get_objective_value()
        best_bound = result.get_objective_bound()

        assignment = []
        station_dict = {}

        # this only collects all tasks in stations
        # however, the tasks are not in the correct order of the sequence
        for s in stations:
            station_dict[s] = {}
            for i in tasks:
                if earliest[i+1] <= s <= latest[i+1]:
                    if result.get_var_solution(t_to_s[s][i]).get_end():
                        station_dict[s][i+1] = result.get_var_solution(t_to_s[s][i]).get_end()
            sorted_station_dict = list(sorted(station_dict[s].items(), key=lambda item: item[1]))
            if sorted_station_dict == []:
                pass
            else:
                assignment.append([ele[0] for ele in sorted_station_dict])

        if result.is_solution_optimal():
            is_optimal = True
            print("optimal cost: {}".format(cost))
        else:
            print("gap: {}".format(result.get_objective_gap()))
            print("best bound: {}".format(result.get_objective_bound()))


    elif result.get_solve_status() == "infeasible":
        is_infeasible = True
        print("The problem is infeasible.")

    return cost, best_bound, is_infeasible, is_optimal, time_used, time_out, assignment, res


        


def compute_m_bounds(number_of_tasks, cycle_time, task_times,
                     forward_min_task, backward_min_task):
    lb = math.ceil(
        (sum(task_times.values()) + 
        sum([min(forward_min_task[i], backward_min_task[i]) for i in range(number_of_tasks)])
        ) / cycle_time
    )
    ub = min(2 * lb, number_of_tasks)
    return lb, ub



def compute_all_predecessors(tasks, predecessors):
    all_predecessors = {}
    for i in tasks:
        predecessor_set = set()
        for j in predecessors[i]:
            predecessor_set.add(j)
            predecessor_set.update(all_predecessors[j])
        all_predecessors[i] = sorted(list(predecessor_set))
    return all_predecessors



def compute_all_followers(tasks, followers):
    all_followers = {}
    for i in reversed(tasks):
        follower_set = set()
        for j in followers[i]:
            follower_set.add(j)
            follower_set.update(all_followers[j])
        all_followers[i] = sorted(list(follower_set))
    return all_followers



def compute_earliest_station(i, cycle_time, task_times, predecessors, 
                             forward_min_task_f, backward_min_task_f):
    return math.ceil(
        (task_times[i] + 
         sum(task_times[j] + 
             min(forward_min_task_f[j-1], 
                 backward_min_task_f[j-1]) 
                 for j in predecessors[i])) / cycle_time
    )



def compute_latest_station(i, cycle_time, task_times, followers, m,
                           forward_min_task_b, backward_min_task_b):
    return (
        m
        + 1
        - math.ceil(
            (task_times[i] + 
             sum(task_times[j] +
                 min(forward_min_task_b[j-1], 
                     backward_min_task_b[j-1]) 
                     for j in followers[i])) / cycle_time
        )
    )



def compute_earliest_latest(number_of_tasks, cycle_time, task_times, predecessors, followers,
                            forward_min_task_f, backward_min_task_f,
                            forward_min_task_b, backward_min_task_b,
                            forward_min_max_f, backward_min_min_f):
    tasks = list(range(1, number_of_tasks + 1))
    m_lb, m_ub = compute_m_bounds(number_of_tasks, cycle_time, task_times,
                                  forward_min_task_f, backward_min_task_f)
    all_predecessors = compute_all_predecessors(tasks, predecessors)
    all_followers = compute_all_followers(tasks, followers)
    earliest = {}
    latest = {}
    for i in tasks:
        earliest[i] = compute_earliest_station(
            i, cycle_time, task_times, all_predecessors,
            forward_min_task_f, backward_min_task_f
        )
        latest[i] = compute_latest_station(
            i, cycle_time, task_times, all_followers, m_ub,
            forward_min_task_b, backward_min_task_b
        )
    return m_lb, m_ub, all_predecessors, all_followers, earliest, latest



# compute station distance between two tasks
def compute_distance(number_of_tasks, task_times, cycle_time, 
                     all_predecessors, all_followers,
                     forward_min_task_f, backward_min_task_f):
    tasks = list(range(1, number_of_tasks+1))
    d = {
        (i, j): math.floor(
            (
                task_times[i]
                + task_times[j]
                + min(forward_min_task_f[i-1], backward_min_task_f[i-1])
                - 1
                + sum(
                    task_times[k]
                    + min(forward_min_task_f[k-1], backward_min_task_f[k-1])
                    for k in all_followers[i]
                    if k in all_predecessors[j]
                )
            )
            / cycle_time
        )
        for j in tasks
        for i in all_predecessors[j]
    }
    return d



def check_sol_sualbp1(sol, number_of_tasks, cycle_time, task_times, 
                      predecessors, followers, forward, backward):
    # 1. check number of tasks
    ct = 0
    for s in sol:
        ct += len(s)
    if ct != number_of_tasks:
        print("Wrong number of tasks in the solution")
        return False

    # 2. check the cycle time at each station
    # 3. check the predecessors for each task
    visited = []
    for s in sol:
        stime = task_times[s[0]]
        for pred in predecessors[s[0]]:
            if pred not in visited:
                print("Processors not sheduled ahead of task {}".format(s[0]))
                return False
        if s[0] not in visited:
            visited.append(s[0])
        if len(s) <= 1:
            pass
        else:
            for i in range(1, len(s)):
                for pred in predecessors[s[i]]:
                    if pred not in visited:
                        print("Processors not sheduled ahead of task {}".format(s[i]))
                        return False
                if s[i] not in visited:
                    visited.append(s[i])
                stime += forward[s[i-1]-1][s[i]-1]
                stime += task_times[s[i]]
            stime += backward[s[len(s)-1]-1][s[0]-1]
        
        # 4. check whether station time is <= cycle time
        if stime > cycle_time:
            print("Wrong station time of station {} with station time {} and cycle time {}".format(s, stime, cycle_time))
            return False
    
    # print("The visited is:", visited)
    assert len(visited) == number_of_tasks, "visited is problematic!"
    
    return True




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str)
    parser.add_argument("--runtime", default=1800, type=int)
    args = parser.parse_args()

    number_of_tasks, cycle_time, task_times, predecessors, followers, forward, backward = \
        utils_sualbsp1.read(args.input)
    
    
    # The indices number_of_tasks and number_of_tasks+1 are dummy for CP
    forward.append([0*i for i in range(number_of_tasks+1)])
    for row in forward:
        row.append(0)
    backward.append([0*i for i in range(number_of_tasks+1)])
    for row in backward:
        row.append(0)

    # for the convenience of element constraint as it supports only the 1D array
    backward_array = []
    for i in range(number_of_tasks+2):
        for j in range(number_of_tasks+2):
            backward_array.append(backward[i][j])
    
    # i -> j
    forward_min_task_f = []
    for i in range(number_of_tasks):
        forward_min_task_ = 1000
        for j in range(number_of_tasks):
            if i!=j:
                forward_min_task_ = min(forward_min_task_, forward[i][j])    # setup to i, instead of from i
        forward_min_task_f.append(forward_min_task_)    
    forward_min_min_f = min(forward_min_task_f)
    forward_min_max_f = max(forward_min_task_f)
    forward_min_task_f.append(0)
    print("forward_min_task_f: ", forward_min_task_f)
    print("forward_min_min_f: ", forward_min_min_f)
    print("forward_min_max_f: ", forward_min_max_f)

    # i -> j
    backward_min_task_f = []
    for i in range(number_of_tasks):
        backward_min_task_ = 1000
        for j in range(number_of_tasks):
            if i!=j:
                backward_min_task_ = min(backward_min_task_, backward[i][j])    # setup to i, instead of from i
        backward_min_task_f.append(backward_min_task_)    
    backward_min_min_f = min(backward_min_task_f)
    backward_min_max_f = max(backward_min_task_f)
    backward_min_task_f.append(0)
    print("backward_min_task_f: ", backward_min_task_f)
    print("backward_min_min_f: ", backward_min_min_f)
    print("backward_min_max_f: ", backward_min_max_f)

    # j -> i
    forward_min_task_b = []
    for i in range(number_of_tasks):
        forward_min_task_ = 1000
        for j in range(number_of_tasks):
            if i!=j:
                forward_min_task_ = min(forward_min_task_, forward[j][i])    # setup to i, instead of from i
        forward_min_task_b.append(forward_min_task_)    
    forward_min_min_b = min(forward_min_task_b)
    forward_min_max_b = max(forward_min_task_b)
    forward_min_task_b.append(0)
    print("forward_min_task_b: ", forward_min_task_b)
    print("forward_min_min_b: ", forward_min_min_b)
    print("forward_min_max_b: ", forward_min_max_b)

    # j -> i
    backward_min_task_b = []
    for i in range(number_of_tasks):
        backward_min_task_ = 1000
        for j in range(number_of_tasks):
            if i!=j:
                backward_min_task_ = min(backward_min_task_, backward[j][i])    # setup to i, instead of from i
        backward_min_task_b.append(backward_min_task_)    
    backward_min_min_b = min(backward_min_task_b)
    backward_min_max_b = max(backward_min_task_b)
    backward_min_task_b.append(0)
    print("backward_min_task_b: ", backward_min_task_b)
    print("backward_min_min_b: ", backward_min_min_b)
    print("backward_min_max_b: ", backward_min_max_b)
    
    # preprocessing
    m_lb, m_ub, all_predecessors, all_followers, earliest, latest = \
        compute_earliest_latest(
            number_of_tasks, cycle_time, task_times, predecessors, followers,
            forward_min_task_f, backward_min_task_f,
            forward_min_task_b, backward_min_task_b,
            forward_min_max_f, backward_min_min_f
        )

    print("number_of_tasks: ", number_of_tasks)
    print("cycle_time: ", cycle_time)
    print("task_times: ", task_times)
    print("predecessors: ", predecessors)
    print("followers: ", followers)
    print("forward: ")
    for row in forward:
        print(row)
    print("backward: ")
    for row in backward:
        print(row)
    for i in range(number_of_tasks+2):
        for j in range(number_of_tasks+2):
            assert backward_array[i*(number_of_tasks+2)+j] == backward[i][j], \
            "backward_array value errors at i={} and j={}!".format(i,j)

    # sys.exit()

    # a lot of preprocessing!
    Ai = utils_sualbsp1.calculate_Ai(number_of_tasks, earliest, latest)

    FiF, PiF, FiB, PiB = utils_sualbsp1.calculate_F_and_P(number_of_tasks, task_times, 
        predecessors, followers,
        all_predecessors, all_followers, 
        forward, backward,
        Ai)
    
    aijF, bijF, aijB, bijB = utils_sualbsp1.calculate_pair_earliest(number_of_tasks, task_times, cycle_time, m_ub,
        all_predecessors, all_followers, forward, backward,
        FiF, PiF, FiB, PiB)

    Bi, Rij = utils_sualbsp1.calculate_B_and_R(number_of_tasks, earliest, latest,
        all_predecessors, all_followers)
    
    delta, deltap = utils_sualbsp1.calculate_delta(number_of_tasks, earliest, latest,
        all_predecessors, all_followers,
        forward, backward,
        Bi, Rij,
        FiF, PiF, FiB, PiB)
    
    dij, AiNew = utils_sualbsp1.calculate_dij(number_of_tasks, task_times, cycle_time,
        earliest, latest,
        all_predecessors, all_followers,
        forward, backward,
        Bi, Rij,
        FiF, PiF, FiB, PiB,
        delta, deltap)
    
    distance = compute_distance(number_of_tasks, task_times, cycle_time, 
        all_predecessors, all_followers,
        forward_min_task_f, backward_min_task_f)

    # cost, best_bound, is_infeasible, is_optimal, time_used, time_out, assignment, res 
    cost, best_bound, is_infeasible, is_optimal, time_used, time_out, assignment, res = solve(
        number_of_tasks,
        cycle_time,
        task_times,
        m_lb, m_ub, all_predecessors, all_followers, earliest, latest,
        forward,
        backward,
        backward_array,
        FiF, PiF, FiB, PiB,
        aijF, bijF, aijB, bijB,
        dij, AiNew,
        distance,
        time_limit=args.runtime,
    )
    if len(assignment) > 0:
        if not check_sol_sualbp1(assignment, number_of_tasks, cycle_time, task_times, 
                                predecessors, followers, forward, backward):
            print("Solution of the instance {} is problematic".format(args.input))
    
    print("*************************************************************")
    print("cost:", cost)
    print("best_bound:", best_bound)
    print("is_infeasible:", is_infeasible)
    print("is_optimal:", is_optimal)
    print("time_used:", time_used)
    print("time_out:", time_out)
    print("assignment:", assignment)
    print("res:", res)

    utils_sualbsp1.write_res("cp", args.input, 
                             cost, best_bound, is_infeasible, is_optimal, time_used, time_out, 
                             assignment, res)

