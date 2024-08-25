import argparse
import math
import time
import sys

import didppy as dp

import utils_sualbsp2


def solve(
    number_of_tasks,
    number_of_stations,
    lb_cycle,
    task_times,
    predecessors,
    followers,
    forward, 
    forward_min_task,
    forward_min_min,
    forward_min_max,
    backward, 
    backward_min_task,
    backward_min_min,
    backward_min_max,
    time_limit=None,
):
    time_start = time.time()

    model = dp.Model()

    task = model.add_object_type(number=number_of_tasks+1)
    uncompleted = model.add_set_var(
        object_type=task, target=list(range(number_of_tasks))
    )
    current_task = model.add_element_var(object_type=task, target=number_of_tasks)
    current_station = model.add_int_var(target=0)
    first_task = model.add_element_var(object_type=task, target=number_of_tasks)
    curr_time = model.add_int_resource_var(target=0, less_is_better=True)
    makespan = model.add_int_resource_var(target=0, less_is_better=True)

    # start from 0, dummy is number_of_tasks+1
    task_time_table = model.add_int_table(
        [task_times[i + 1] for i in range(number_of_tasks)]
    )

    # start from 0, dummy is number_of_tasks+1
    forward_min_task_table = model.add_int_table(
        [forward_min_task[i] for i in range(number_of_tasks+1)]
    )
    backward_min_task_table = model.add_int_table(
        [backward_min_task[i] for i in range(number_of_tasks+1)]
    )

    # start from 0, dummy is number_of_tasks+1
    forward_table = model.add_int_table(forward)
    backward_table = model.add_int_table(backward)

    predecessors_table = model.add_set_table(
        [[j - 1 for j in predecessors[i + 1]] for i in range(number_of_tasks)],
        object_type=task,
    )

    model.add_base_case([
        uncompleted.is_empty(), 
        first_task==number_of_tasks, 
        # current_station==number_of_stations,
    ])

    name_to_task = {}

    # assign job i+1 to a new station
    # update cost by 1
    for i in range(number_of_tasks):
        name = "assign_first {}".format(i+1)
        name_to_task[name] = i + 1
        t = dp.Transition(
            name = name,
            cost = dp.max(dp.IntExpr.state_cost(), task_time_table[i]),
            effects = [
                (uncompleted, uncompleted.remove(i)),
                (curr_time, task_time_table[i]),
                (current_task, i),
                (first_task, i),
                (current_station, current_station+1),
                (makespan, dp.max(makespan, task_time_table[i])),
            ],
            preconditions = [
                uncompleted.contains(i),
                (uncompleted & predecessors_table[i]).is_empty(),
                first_task == number_of_tasks,
                current_station < number_of_stations,
            ]
        )
        model.add_transition(t)

    # schedule i+1 as the next task of the current station
    for i in range(number_of_tasks):
        name = "schedule {}".format(i+1)
        name_to_task[name] = i + 1
        t = dp.Transition(
            name=name,
            cost=dp.max(dp.IntExpr.state_cost(), curr_time + task_time_table[i] + forward_table[current_task, i]),
            effects=[
                (uncompleted, uncompleted.remove(i)),
                (curr_time, curr_time + task_time_table[i] + forward_table[current_task, i]),
                (current_task, i),
                (makespan, dp.max(makespan, curr_time + task_time_table[i] + forward_table[current_task, i])),
            ],
            preconditions=[
                first_task < number_of_tasks,
                uncompleted.contains(i),
                (uncompleted & predecessors_table[i]).is_empty(),
            ],
        )
        model.add_transition(t)

    t = dp.Transition(
        name="close the current station",
        cost=dp.max(dp.IntExpr.state_cost(), curr_time + backward_table[current_task, first_task]),
        effects=[
            (curr_time, curr_time + backward_table[current_task, first_task]),
            (makespan, dp.max(makespan, curr_time + backward_table[current_task, first_task])),
            (current_task, number_of_tasks),
            (first_task, number_of_tasks),
        ],
        preconditions=[
            # current_station <= number_of_stations,
            first_task < number_of_tasks,
        ],
    )
    model.add_transition(t)

    model.add_dual_bound(
        math.ceil(
            (task_time_table[uncompleted] + 
                forward_min_task_table[uncompleted] + 
                curr_time - 
                (uncompleted.is_empty()).if_then_else(0,forward_min_task_table.max(uncompleted))*(number_of_stations-current_station) +
                (uncompleted.is_empty()).if_then_else(0, backward_min_task_table.min(uncompleted))*(number_of_stations-current_station) +
                backward_min_task_table[first_task] )
            / dp.min(number_of_stations, number_of_stations-current_station+1)
        ) - makespan
    )

    model.add_dual_bound(
        math.ceil(
            (task_time_table[uncompleted] + curr_time + backward_min_task_table[first_task])
            / dp.min(number_of_stations, number_of_stations-current_station+1)
        ) - makespan
    )

    model.add_dual_bound(0)


    res = []

    solver = dp.CABS(model, time_limit=time_limit, quiet=False)
    # solver = dp.CAASDy(model, time_limit=time_limit, quiet=False)
    # solution = solver.search()
    terminated = False
    time_end = time.time()
    while not terminated and time_end-time_start < time_limit:
        solution, terminated = solver.search_next()
        time_end = time.time()
        # print('The runtime of DIDP is: ', time_end-time_start, 's')
        res.append((time_end-time_start, solution.cost, solution.best_bound))

    if solution.is_infeasible:
        print("Problem is infeasible")
    elif solution.is_optimal:
        print("Optimally solved")
    elif solution.time_out:
        print("Time out")

    assignment = [[]]

    for t in solution.transitions:
        print(t.name)

        if t.name == "close the current station":
            assignment.append([])
        else:
            assignment[-1].append(name_to_task[t.name])

    print("expanded: {}".format(solution.expanded))
    print("cost: {}".format(solution.cost))

    return solution.cost, solution.best_bound, \
           solution.is_infeasible, solution.is_optimal, \
           solution.time, solution.time_out, \
           assignment, res



def compute_m_bounds(number_of_tasks, cycle_time, task_times,
                     forward_min_task, backward_min_task):
    lb = math.ceil(
        (sum(task_times.values()) + 
        sum([min(forward_min_task[i], backward_min_task[i]) for i in range(number_of_tasks)])
        ) / cycle_time
    )
    ub = min(2 * lb, number_of_tasks)
    return lb, ub


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
   
    forward_min_task = []
    for i in range(number_of_tasks):
        forward_min_task_ = 1000
        for j in range(number_of_tasks):
            if i!=j:
                forward_min_task_ = min(forward_min_task_, forward[j][i])    # setup to i, instead of from i
        forward_min_task.append(forward_min_task_)    
    forward_min_min = min(forward_min_task)
    forward_min_max = max(forward_min_task)
    forward_min_task.append(0)
    print("forward_min_task: ", forward_min_task)
    print("forward_min_min: ", forward_min_min)
    print("forward_min_max: ", forward_min_max)

    backward_min_task = []
    for i in range(number_of_tasks):
        backward_min_task_ = 1000
        for j in range(number_of_tasks):
            if i!=j:
                backward_min_task_ = min(backward_min_task_, backward[j][i])    # setup to i, instead of from i
        backward_min_task.append(backward_min_task_)    
    backward_min_min = min(backward_min_task)
    backward_min_max = max(backward_min_task)
    backward_min_task.append(0)
    print("backward_min_task: ", backward_min_task)
    print("backward_min_min: ", backward_min_min)
    print("backward_min_max: ", backward_min_max)

    # m_lb, m_ub = compute_m_bounds(number_of_tasks, cycle_time, task_times, 
    #                               forward_min_task, backward_min_task)

    lb_cycle, forward_min_task_sorted, backward_min_task_sorted = \
        compute_lb_cycle_time(number_of_tasks, number_of_stations, task_times, forward_min_task, backward_min_task)
    print("The lower bound of cycle time is: ", lb_cycle)


    print("number_of_tasks: ", number_of_tasks)
    print("number_of_stations:", number_of_stations)
    print("lb_cycle: ", lb_cycle)
    print("task_times: ", task_times)
    print("predecessors: ", predecessors)
    print("followers: ", followers)

    
    cost, best_bound, is_infeasible, is_optimal, time_used, time_out, assignment, res = solve(
        number_of_tasks,
        number_of_stations,
        lb_cycle,
        task_times,
        predecessors,
        followers,
        forward, 
        forward_min_task,
        forward_min_min,
        forward_min_max,
        backward, 
        backward_min_task,
        backward_min_min,
        backward_min_max,
        time_limit=args.runtime,
    )
    
    print("*************************************************************")
    print("cost:", cost)
    print("best_bound:", best_bound)
    print("is_infeasible:", is_infeasible)
    print("is_optimal:", is_optimal)
    print("time_used:", time_used)
    print("time_out:", time_out)
    print("assignment:", assignment)
    print("res:", res)



    utils_sualbsp2.write_res("cabs", args.input, 
                             cost, best_bound, is_infeasible, is_optimal, time_used, time_out, 
                             assignment, res)

