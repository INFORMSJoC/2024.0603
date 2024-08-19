import argparse
import math
import time
import sys

import didppy as dp

import utils_sualbsp1


def solve(
    number_of_tasks,
    m_lb, m_ub,
    cycle_time,
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
    idle_time = model.add_int_resource_var(target=0, less_is_better=False)
    current_task = model.add_element_var(object_type=task, target=number_of_tasks)
    current_station = model.add_int_resource_var(target=0, less_is_better=True)
    first_task = model.add_element_var(object_type=task, target=number_of_tasks)

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

    lb2_weight1 = model.add_int_table(
        [1 if task_times[i + 1] > cycle_time / 2 else 0 for i in range(number_of_tasks)]
    )
    lb2_weight2 = model.add_float_table(
        [
            0.5 if task_times[i + 1] == cycle_time / 2 else 0
            for i in range(number_of_tasks)
        ]
    )
    lb3_weight = model.add_float_table(
        [
            1.0
            if task_times[i + 1] > cycle_time * 2 / 3
            else 2 / 3 // 0.001 / 1000
            if task_times[i + 1] == cycle_time * 2 / 3
            else 0.5
            if task_times[i + 1] > cycle_time / 3
            else 1 / 3 // 0.001 / 1000
            if task_times[i + 1] == cycle_time / 3
            else 0.0
            for i in range(number_of_tasks)
        ]
    )

    model.add_base_case([
        uncompleted.is_empty(), 
        first_task==number_of_tasks, 
    ])

    name_to_task = {}

    # assign job i+1 to a new station
    # update cost by 1
    for i in range(number_of_tasks):
        name = "assign_first {}".format(i+1)
        name_to_task[name] = i + 1
        t = dp.Transition(
            name = name,
            cost = dp.IntExpr.state_cost() + 1,
            effects = [
                (uncompleted, uncompleted.remove(i)),
                (idle_time, cycle_time - task_time_table[i]),
                (current_task, i),
                (first_task, i),
                (current_station, current_station+1),
            ],
            preconditions = [
                uncompleted.contains(i),
                (uncompleted & predecessors_table[i]).is_empty(),
                first_task == number_of_tasks,
            ]
        )
        model.add_transition(t)

    # schedule i+1 as the next task of the current station
    for i in range(number_of_tasks):
        name = "schedule {}".format(i+1)
        name_to_task[name] = i + 1
        t = dp.Transition(
            name=name,
            cost=dp.IntExpr.state_cost(),
            effects=[
                (uncompleted, uncompleted.remove(i)),
                (idle_time, idle_time - task_time_table[i] - forward_table[current_task, i]),
                (current_task, i),
            ],
            preconditions=[
                first_task < number_of_tasks,
                uncompleted.contains(i),
                task_time_table[i] + forward_table[current_task, i] <= idle_time,
                (uncompleted & predecessors_table[i]).is_empty(),
            ],
        )
        model.add_transition(t)

    t = dp.Transition(
        name="close the current station",
        cost=dp.IntExpr.state_cost(),
        effects=[
            (idle_time, 0),
            (current_task, number_of_tasks),
            (first_task, number_of_tasks),
        ],
        preconditions=[
            ~uncompleted.contains(i)
            | (task_time_table[i] + forward_table[current_task, i] + backward_table[i, first_task] > idle_time)
            | ((uncompleted & predecessors_table[i]).len() > 0)
            for i in range(number_of_tasks)
        ] + [
            first_task < number_of_tasks,
            backward_table[current_task, first_task] <= idle_time,
        ],
    )
    model.add_transition(t)

    model.add_dual_bound(
        math.ceil((task_time_table[uncompleted] - idle_time) / cycle_time)
    )

    model.add_dual_bound(
        math.ceil(
            (forward_min_task_table[uncompleted] + 
             task_time_table[uncompleted] + 
             (((m_lb - current_station) <= 0) | uncompleted.is_empty()).if_then_else(0, 
                (m_lb - current_station) * backward_min_task_table.min(uncompleted)
             ) +
             backward_min_task_table[first_task] - 
             (m_ub - current_station) * (uncompleted.is_empty()).if_then_else(0,
                forward_min_task_table.max(uncompleted)
             ) - 
             idle_time) / cycle_time
        )
    )

    model.add_dual_bound(
        lb2_weight1[uncompleted]
        + math.ceil(lb2_weight2[uncompleted])
        - (idle_time >= cycle_time / 2).if_then_else(1, 0)
    )
    model.add_dual_bound(
        math.ceil(lb3_weight[uncompleted])
        - (idle_time >= cycle_time / 3).if_then_else(1, 0)
    )

    model.add_dual_bound(0)

    res = []

    
    solver = dp.CABS(model, time_limit=time_limit, quiet=False)
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



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str)
    parser.add_argument("--runtime", default=1800, type=int)
    args = parser.parse_args()

    number_of_tasks, cycle_time, task_times, predecessors, followers, forward, backward = \
        utils_sualbsp1.read(args.input)
        
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

    m_lb, m_ub = compute_m_bounds(number_of_tasks, cycle_time, task_times, 
                                  forward_min_task, backward_min_task)


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

    
    cost, best_bound, is_infeasible, is_optimal, time_used, time_out, assignment, res = solve(
        number_of_tasks,
        m_lb, m_ub,
        cycle_time,
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

    utils_sualbsp1.write_res("cabs", args.input, 
                             cost, best_bound, is_infeasible, is_optimal, time_used, time_out, 
                             assignment, res)
