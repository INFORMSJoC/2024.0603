import argparse
import math
import time

import gurobipy as gp

import utils_sualbsp2


def solve(
    number_of_tasks,
    number_of_stations,
    lb_cycle_new,
    ub_cycle_new,
    task_times,
    m_lb, m_ub, all_predecessors, all_followers, earliest, latest,
    forward,
    backward,
    FiF, PiF, FiB, PiB,
    aijF, bijF, aijB, bijB,
    dij, AiNew,
    time_limit=None,
    threads=1,
    verbose=False,
    history=None,
):
    tasks = list(range(1, number_of_tasks + 1))
    stations = list(range(1, number_of_stations + 1))

    x_indices = []
    for i in tasks:
        x_indices += [(s, i) for s in range(earliest[i], latest[i] + 1)]
      
    g_indices = []
    for s in stations:
        for i in range(1, number_of_tasks+1):
            if s in range(earliest[i], latest[i]+1):
                for j in range(1, number_of_tasks+1):
                    if s in range(earliest[j], latest[j]+1):
                        if j in FiF[i]:
                            if s in range(aijF[(i,j)], bijF[(i,j)]+1):
                                if i!=j:
                                    g_indices.append((i,j,s))

    h_indices = []
    for s in stations:
        for i in range(1, number_of_tasks+1):
            if s in range(earliest[i], latest[i]+1):
                for j in range(1, number_of_tasks+1):
                    if s in range(earliest[j], latest[j]+1):
                        if j in FiB[i]:
                            if s in range(aijB[(i,j)], bijB[(i,j)]+1):
                                h_indices.append((i,j,s))
            if s in range(earliest[i], latest[i]+1) and (i,i,s) not in h_indices:
                h_indices.append((i,i,s))
            
    model = gp.Model()
    model.setParam("Threads", threads)
    if time_limit is not None:
        model.setParam("TimeLimit", time_limit)
    # if not verbose:
    #     model.setParam("OutputFlag", 0)


    ####################### Decision Variables ###############################
    # assignment
    x = model.addVars(x_indices, vtype=gp.GRB.BINARY)

    # sequence of tasks on stations
    g = model.addVars(g_indices, vtype=gp.GRB.BINARY)

    # last-first tasks on stations
    h = model.addVars(h_indices, vtype=gp.GRB.BINARY)

    # auxiliary for sub-tour elimination
    r = model.addVars(tasks, vtype=gp.GRB.INTEGER, lb=0)

    # the station of tasks
    z = model.addVars(tasks, vtype=gp.GRB.INTEGER, lb=1)

    makespan = model.addVar(vtype=gp.GRB.CONTINUOUS, obj=1, 
                            lb=lb_cycle_new, ub=ub_cycle_new)

    # ORS (2)
    model.addConstrs(
        gp.quicksum(x[s, i]   
                    for s in stations 
                    if earliest[i] <= s <= latest[i]) 
        == 1
        for i in tasks
    )

    # IJOC (25)
    model.addConstrs(
        makespan + 
        ub_cycle_new *
        (gp.quicksum(s * x[s,j] for s in range(earliest[j], latest[j]+1)) - 
         gp.quicksum(s * x[s,i] for s in range(earliest[i], latest[i]+1)))
        >= dij[(i,j)]
        for i in tasks
        for j in FiF[i]
        if j not in AiNew[i]
        if (i,j) in dij.keys()
    )

    # IJOC (26)
    model.addConstrs(
        gp.quicksum(x[s, i]   
                    for i in tasks 
                    if earliest[i] <= s <= latest[i]) 
        >= 1
        for s in stations
    )

    # ORS (3)
    model.addConstrs(
        gp.quicksum(x[s, i] * s 
                    for s in stations 
                    if earliest[i] <= s <= latest[i]) 
        == z[i]
        for i in tasks
    )

    # ORS (39)
    model.addConstrs(
        gp.quicksum(g[i,j,s]
                    for j in tasks 
                    if earliest[j] <= s <= latest[j]
                    if (i,j,s) in g_indices)
        +
        gp.quicksum(h[i,j,s]
                    for j in tasks 
                    if earliest[j] <= s <= latest[j]
                    if (i,j,s) in h_indices)
        == x[s, i]
        for i in tasks
        for s in stations
        if earliest[i] <= s <= latest[i]
    )

    # ORS (40)
    model.addConstrs(
        gp.quicksum(g[i,j,s]
                    for i in tasks 
                    if earliest[i] <= s <= latest[i]
                    if (i,j,s) in g_indices)
        +
        gp.quicksum(h[i,j,s]
                    for i in tasks 
                    if earliest[i] <= s <= latest[i]
                    if (i,j,s) in h_indices)
        == x[s, j]
        for j in tasks
        for s in stations
        if earliest[j] <= s <= latest[j]
    )

    # ORS (41)
    # it is possible that h[(i,i,s)]==1
    # it means that i is the first task of station s
    model.addConstrs(
        gp.quicksum(h[i,j,s] 
                    for i in tasks 
                    for j in tasks 
                    if (i,j,s) in h_indices) == 1
        for s in stations
    )

    # ORS (43)
    model.addConstrs(
        r[i] + 1 + 
        number_of_tasks * (gp.quicksum(g[i,j,s] 
                            for s in stations 
                            if (i,j,s) in g_indices) - 1) 
        <= r[j]
        for i in tasks 
        for j in FiF[i]
    )

    # ORS (44)
    model.addConstrs(
        r[i] + 1 <= r[j]
        for i in tasks for j in all_followers[i]
    )

    # ORS (45)
    model.addConstrs(
        z[i] <= z[j]
        for i in tasks for j in all_followers[i]
    )

    # ORS (46)
    model.addConstrs(
        gp.quicksum(g[i,j,s] * forward[i-1][j-1]
                    for i in tasks
                    for j in tasks 
                    if earliest[i] <= s <= latest[i]
                    if earliest[j] <= s <= latest[j]
                    if (i,j,s) in g_indices) 
        +
        gp.quicksum(h[i,j,s] * backward[i-1][j-1]
                    for i in tasks
                    for j in tasks 
                    if earliest[i] <= s <= latest[i]
                    if earliest[j] <= s <= latest[j]
                    if (i,j,s) in h_indices) 
        +
        gp.quicksum(x[s,i] * task_times[i]
                    for i in tasks
                    if earliest[i] <= s <= latest[i]) 
        <= makespan
        for s in stations
    )
    
    # ORS (48)
    model.addConstrs(
        gp.quicksum(x[s, i] 
                    for i in tasks
                    if i!=j
                    if earliest[i] <= s <= latest[i]) 
        <= (number_of_tasks - number_of_stations + 1) * (1 - h[j,j,s])
        for j in tasks
        for s in stations
        if earliest[j] <= s <= latest[j]
        if (j,j,s) in h_indices
    )

    # ORS (51)
    model.addConstrs(
        r[i] <= number_of_tasks - len(all_followers[i])
        for i in tasks
    )
    model.addConstrs(
        r[i] >= 1 + len(all_predecessors[i])
        for i in tasks
    )

    time_start = time.time()

    def get_callback(res):
        def dump_solution(model, where):
            if where == gp.GRB.Callback.MIPSOL:
                res.append(
                        (time.time() - time_start,
                        model.cbGet(gp.GRB.Callback.MIPSOL_OBJ),
                        model.cbGet(gp.GRB.Callback.MIPSOL_OBJBND))
                )

                if time.time() - time_start > time_limit:
                    model.terminate()

        return dump_solution

    is_infeasible = False
    is_optimal = False
    cost = -1
    best_bound = -1
    assignment = []
    res = []
    callback = get_callback(res)
    model.optimize(callback)

    time_used = time.time() - time_start
    time_out = time.time() - time_start

    status = model.getAttr("Status")
    sol_count = model.getAttr("SolCount")
    if status == gp.GRB.INFEASIBLE:
        print("infeasible")
        is_infeasible = True
        return cost, best_bound, is_infeasible, is_optimal, time_used, time_out, assignment, res
    elif sol_count > 0:
        solution = []
        for s in stations:
            tasks_in_station = []
            for i in tasks:
                if earliest[i] <= s <= latest[i] and x[s, i].X > 0.5:
                    tasks_in_station.append(i)
            if len(tasks_in_station) > 0:
                solution.append(tasks_in_station.copy())
        cost = round(model.objVal)
        print("solution: ")
        for s in solution:
            print(s)
        print("cost: {}".format(cost))

        # solution of the rank
        sol_rank = {}
        for i in tasks:
            sol_rank[i] = int(r[i].X)
        
        # sort tasks in stations
        solution_sorted = []
        sol_rank_station_dict = {}
        for s in range(len(solution)):
            sol_rank_station_dict[s] = {}
            for i in solution[s]:
                sol_rank_station_dict[s][i] = sol_rank[i]
            sorted_station_dict = list(sorted(sol_rank_station_dict[s].items(), 
                                              key=lambda item: item[1]))
            if sorted_station_dict == []:
                pass
            else:
                solution_sorted.append([ele[0] for ele in sorted_station_dict])
        
        print("solution_sorted: ")
        for s in solution_sorted:
            print(s)

        best_bound = model.ObjBound

        if status == gp.GRB.OPTIMAL:
            print("optimal cost: {}".format(cost))
            is_optimal = True
        else:
            print("gap: {}".format(model.getAttr("MIPGap")))
            print("best bound: {}".format(model.getAttr("ObjBound")))
        return cost, best_bound, is_infeasible, is_optimal, time_used, time_out, solution_sorted, res
    else:
        return cost, best_bound, is_infeasible, is_optimal, time_used, time_out, assignment, res
    return cost, best_bound, is_infeasible, is_optimal, time_used, time_out, assignment, res
        



def compute_m_bounds(number_of_tasks, ub_cycle, task_times,
                     forward_min_task, backward_min_task):
    lb = math.ceil(
        (sum(task_times.values()) + 
        sum([min(forward_min_task[i], backward_min_task[i]) for i in range(number_of_tasks)])
        ) / ub_cycle
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



def compute_earliest_station(i, ub_cycle, task_times, predecessors, 
                             forward_min_task_f, backward_min_task_f):
    return math.ceil(
        (task_times[i] + 
         sum(task_times[j] + 
             min(forward_min_task_f[j-1], 
                 backward_min_task_f[j-1]) 
                 for j in predecessors[i])) / ub_cycle
    )



def compute_latest_station(i, ub_cycle, task_times, followers, m,
                           forward_min_task_b, backward_min_task_b):
    return (
        m
        + 1
        - math.ceil(
            (task_times[i] + 
             sum(task_times[j] +
                 min(forward_min_task_b[j-1], 
                     backward_min_task_b[j-1]) 
                     for j in followers[i])) / ub_cycle
        )
    )



def compute_earliest_latest(number_of_tasks, number_of_stations, 
                            ub_cycle, task_times, predecessors, followers,
                            forward_min_task_f, backward_min_task_f,
                            forward_min_task_b, backward_min_task_b,
                            forward_min_max_f, backward_min_min_f):
    tasks = list(range(1, number_of_tasks + 1))
    m_lb, m_ub = compute_m_bounds(number_of_tasks, ub_cycle, task_times,
                                  forward_min_task_f, backward_min_task_f)
    all_predecessors = compute_all_predecessors(tasks, predecessors)
    all_followers = compute_all_followers(tasks, followers)
    earliest = {}
    latest = {}
    for i in tasks:
        earliest[i] = compute_earliest_station(
            i, ub_cycle, task_times, all_predecessors,
            forward_min_task_f, backward_min_task_f
        )
        latest[i] = compute_latest_station(
            i, ub_cycle, task_times, all_followers, number_of_stations,
            forward_min_task_b, backward_min_task_b
        )
    return m_lb, m_ub, all_predecessors, all_followers, earliest, latest



def compute_lb_cycle_time(number_of_tasks, number_of_stations, task_times, forward_min_task, backward_min_task):
    '''
    Compute a lower bound of cycle time
    '''
    forward_min_task_sorted = sorted(forward_min_task)
    total = 0
    for task in task_times.values():
        total += task
    for i in range(number_of_tasks - number_of_stations):
        total += forward_min_task_sorted[i]

    backward_min_task_sorted = sorted(forward_min_task)
    for i in range(number_of_stations):
        total += backward_min_task_sorted[i]
    
    return int(math.ceil(total / number_of_stations)), forward_min_task_sorted, backward_min_task_sorted




def compute_ub_cycle_time(number_of_tasks, number_of_stations, task_times, forward, 
                          forward_max_task, backward_max_task):
    '''
    Compute an upper bound of cycle time
    '''
    option1 = 0
    for i in range(number_of_tasks):
        if task_times[i+1] + forward[i][i] > option1:
            option1 = task_times[i+1] + forward[i][i]

    option2 = 0
    for task in task_times.values():
        option2 += task
    forward_max_task_sorted = sorted(forward_max_task, reverse=True)
    for i in range(number_of_tasks - number_of_stations):
        option2 += forward_max_task_sorted[i]
    backward_max_task_sorted = sorted(backward_max_task, reverse=True)
    for i in range(number_of_stations):
        option2 += backward_max_task_sorted[i]
    
    if number_of_stations // 2 == 0:
        option2 = math.floor(
            2*(option2-1) / number_of_stations
        )
    else:
        option2 = math.floor(
            2*(option2) / (number_of_stations+1)
        )

    return max(option1, option2)



##### Need a solution checker!
def check_sol_sualbp2(sol, number_of_tasks, number_of_stations, 
                      cycle_time, task_times, predecessors, followers, forward, backward):
    # 0. check number of stations
    if len(sol) != number_of_stations:
        print("Wrong number of stations in the solution")
        return False
    
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
        utils_sualbsp2.read(args.input)
    
    mlist = utils_sualbsp2.read_xls()

    number_of_stations = mlist[args.input.split('/')[-1].lower()]
    print('number_of_stations:', number_of_stations)
    
    # i -> j
    forward_min_task_f = []
    forward_max_task = []
    for i in range(number_of_tasks):
        forward_min_task_ = 1000
        forward_max_task_ = 0
        for j in range(number_of_tasks):
            if i!=j:
                forward_min_task_ = min(forward_min_task_, forward[i][j])    # setup to i, instead of from i
                forward_max_task_ = max(forward_max_task_, forward[i][j])
        forward_min_task_f.append(forward_min_task_) 
        forward_max_task.append(forward_max_task_)
    forward_min_min_f = min(forward_min_task_f)
    forward_min_max_f = max(forward_min_task_f)
    forward_min_task_f.append(0)
    forward_max_task.append(0)
    print("forward_min_task_f: ", forward_min_task_f)
    print("forward_min_min_f: ", forward_min_min_f)
    print("forward_min_max_f: ", forward_min_max_f)

    # i -> j
    backward_min_task_f = []
    backward_max_task = []
    for i in range(number_of_tasks):
        backward_min_task_ = 1000
        backward_max_task_ = 0
        for j in range(number_of_tasks):
            if i!=j:
                backward_min_task_ = min(backward_min_task_, backward[i][j])    # setup to i, instead of from i
                backward_max_task_ = max(backward_max_task_, backward[i][j])
        backward_min_task_f.append(backward_min_task_)    
        backward_max_task.append(backward_max_task_)
    backward_min_min_f = min(backward_min_task_f)
    backward_min_max_f = max(backward_min_task_f)
    backward_min_task_f.append(0)
    backward_max_task.append(0)
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
    
    lb_cycle, forward_min_task_sorted, backward_min_task_sorted = \
        compute_lb_cycle_time(number_of_tasks, number_of_stations, task_times, 
                              forward_min_task_b, backward_min_task_b)
    print("The lower bound of cycle time is: ", lb_cycle)

    # Use an upper bound of the cycle time to compute the ealiest and latest
    ub_cycle = \
    compute_ub_cycle_time(number_of_tasks, number_of_stations, task_times, forward, 
                          forward_max_task, backward_max_task)
    
    m_lb, m_ub, all_predecessors, all_followers, earliest, latest = \
    compute_earliest_latest(
        number_of_tasks, number_of_stations,
        ub_cycle, task_times, predecessors, followers,
        forward_min_task_f, backward_min_task_f,
        forward_min_task_b, backward_min_task_b,
        forward_min_max_f, backward_min_min_f
    )
    

    # a lot of preprocessing!
    Ai = utils_sualbsp2.calculate_Ai(number_of_tasks, earliest, latest)

    FiF, PiF, FiB, PiB = utils_sualbsp2.calculate_F_and_P(number_of_tasks, task_times, 
        predecessors, followers,
        all_predecessors, all_followers, 
        forward, backward,
        Ai)
    
    aijF, bijF, aijB, bijB = utils_sualbsp2.calculate_pair_earliest(number_of_tasks, task_times, ub_cycle, number_of_stations,
        all_predecessors, all_followers, forward, backward,
        FiF, PiF, FiB, PiB)

    Bi, Rij = utils_sualbsp2.calculate_B_and_R(number_of_tasks, earliest, latest,
        all_predecessors, all_followers)
    
    delta, deltap = utils_sualbsp2.calculate_delta(number_of_tasks, earliest, latest,
        all_predecessors, all_followers,
        forward, backward,
        Bi, Rij,
        FiF, PiF, FiB, PiB)
    
    dij, AiNew = utils_sualbsp2.calculate_dij(number_of_tasks, task_times, ub_cycle,
        earliest, latest,
        all_predecessors, all_followers,
        forward, backward,
        Bi, Rij,
        FiF, PiF, FiB, PiB,
        delta, deltap)

    gamma_i1, gamma_i2, sdF_sorted, sdB_sorted = \
        utils_sualbsp2.calculate_gamma(number_of_tasks, forward, backward,
            FiF, PiF, FiB, PiB)
    
    lb_cycle_new = utils_sualbsp2.calculate_lb_cycle_new(number_of_tasks,   
        number_of_stations, task_times, forward, backward, gamma_i1, gamma_i2, sdF_sorted, sdB_sorted)
    
    ub_cycle_new = utils_sualbsp2.compute_ub_cycle_time_new(
        number_of_tasks, number_of_stations, task_times,
        all_predecessors, all_followers,
        forward, backward, gamma_i1, gamma_i2,
        FiF, PiF, FiB, PiB
    )
    
    earliestn, latestn = \
        utils_sualbsp2.calculate_new_earliest_and_latest(
            number_of_tasks, number_of_stations, task_times, ub_cycle_new,
            all_predecessors, all_followers,
            forward, backward, gamma_i1, gamma_i2)
    
    for i in range(number_of_tasks):
        earliest[i+1] = max(earliest[i+1], earliestn[i+1])
        latest[i+1] = min(latest[i+1], latestn[i+1])
    
    aijF, bijF, aijB, bijB = utils_sualbsp2.calculate_pair_earliest(number_of_tasks, task_times, ub_cycle_new, number_of_stations,
        all_predecessors, all_followers, forward, backward,
        FiF, PiF, FiB, PiB)
    
    dij, AiNew = utils_sualbsp2.calculate_dij(number_of_tasks, task_times, 
        ub_cycle_new,
        earliest, latest,
        all_predecessors, all_followers,
        forward, backward,
        Bi, Rij,
        FiF, PiF, FiB, PiB,
        delta, deltap)


    print("number_of_tasks: ", number_of_tasks)
    print("lb_cycle: ", lb_cycle)
    print("ub_cycle: ", ub_cycle)
    print("lb_cycle_new: ", lb_cycle_new)
    print("ub_cycle_new: ", ub_cycle_new)
    print("task_times: ", task_times)
    print("predecessors: ", predecessors)
    print("followers: ", followers)

    cost, best_bound, is_infeasible, is_optimal, time_used, time_out, assignment, res = solve(
        number_of_tasks,
        number_of_stations,
        lb_cycle_new,
        ub_cycle_new,
        task_times,
        m_lb, m_ub, all_predecessors, all_followers, earliest, latest,
        forward,
        backward,
        FiF, PiF, FiB, PiB,
        aijF, bijF, aijB, bijB,
        dij, AiNew,
        time_limit=args.runtime,
    )
    
    if not check_sol_sualbp2(assignment, number_of_tasks, number_of_stations, cost, task_times, predecessors, followers, forward, backward):
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

    utils_sualbsp2.write_res("mip", args.input, 
                             cost, best_bound, is_infeasible, is_optimal, time_used, time_out, 
                             assignment, res)

