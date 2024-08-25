import os
import pandas as pd
import math


def read(filename):
    print("filename:", filename)

    with open(filename) as f:
        line = f.readline().strip()
        while line != "<end>":
            if line == "<number of tasks>":
                number_of_tasks = int(f.readline().strip())
                # initialize forward and backward setup time
                forward = {
                    (i,j): 0 
                    for i in range(1, number_of_tasks + 1) 
                    for j in range(1, number_of_tasks + 1)
                }
                backward = {
                    (i,j): 0 
                    for i in range(1, number_of_tasks + 1) 
                    for j in range(1, number_of_tasks + 1)
                }
            if line == "<cycle time>":
                cycle_time = int(f.readline().strip())
            if line == "<task times>":
                task_times = {}
                for _ in range(number_of_tasks):
                    pair = f.readline().strip().split()
                    task_times[int(pair[0])] = int(pair[1])
            if line == "<precedence relations>":
                predecessors = {i: [] for i in range(1, number_of_tasks + 1)}
                followers = {i: [] for i in range(1, number_of_tasks + 1)}
                relation = f.readline().strip()
                while len(relation) > 0:
                    relation = relation.split(",")
                    predecessor = int(relation[0])
                    follower = int(relation[1])
                    predecessors[follower].append(predecessor)
                    followers[predecessor].append(follower)
                    relation = f.readline().strip()
            if line == "<setup times forward>":
                # print("forward:", forward)
                setup = f.readline().strip()
                while len(setup) > 0:
                    setup1 = setup.split(",")
                    i = int(setup1[0])
                    setup2 = setup1[1].split(":")
                    j = int(setup2[0])
                    v = int(setup2[1])
                    forward[(i,j)] = v
                    setup = f.readline().strip()
            if line == "<setup times backward>":
                setup = f.readline().strip()
                while len(setup) > 0:
                    setup1 = setup.split(",")
                    i = int(setup1[0])
                    setup2 = setup1[1].split(":")
                    j = int(setup2[0])
                    v = int(setup2[1])
                    backward[(i,j)] = v
                    setup = f.readline().strip()
            
            line = f.readline().strip()

    forward_matrix = []
    backward_matrix = []
    for i in range(1, number_of_tasks + 1):
        forward_matrix_ = []
        backward_matrix_ = []
        for j in range(1, number_of_tasks + 1):
            forward_matrix_.append(forward[(i,j)])
            backward_matrix_.append(backward[(i,j)])
        forward_matrix_.append(0)
        backward_matrix_.append(0)
        forward_matrix.append(forward_matrix_)
        backward_matrix.append(backward_matrix_)
    forward_matrix.append([0*i for i in range(number_of_tasks+1)])
    backward_matrix.append([0*i for i in range(number_of_tasks+1)])

    return number_of_tasks, cycle_time, task_times, predecessors, followers, forward_matrix, backward_matrix






def write_res(header, filename, cost, best_bound, is_infeasible, is_optimal, time_used, time_out, assignment, res):
    
    sep = filename.split('/')
    print(sep)

    folder = "res2/" + header + "/" + sep[2] + "/" + sep[3] + "/"
    os.makedirs(folder, exist_ok=True)

    new_filename = "res2/" + header + "/" + sep[2] + "/" + sep[3] + "/" + sep[4][:-3] + "stats"
    print(new_filename)

    f = open(new_filename, 'w')

    f.write(str(cost) + " ")
    f.write(str(best_bound) + " ")
    f.write(str(is_infeasible) + " ")
    f.write(str(is_optimal) + " ")
    f.write(str(time_used) + " ")
    f.write(str(time_out) + "\n")

    for row in assignment:
        for i in row:
            f.write(str(i) + " ")
        f.write("\n")
    
    for row in res:
        f.write(str(row[0]) + " ")
        f.write(str(row[1]) + " ")
        f.write(str(row[2]) + " ")
        f.write("\n")




def write_time(header, filename, time_used_list):
    sep = filename.split('/')
    print(sep)

    folder = "times/" + header + "/" + sep[2] + "/" + sep[3] + "/"
    os.makedirs(folder, exist_ok=True)

    new_filename = "times/" + header + "/" + sep[2] + "/" + sep[3] + "/" + sep[4][:-3] + "times"
    print(new_filename)

    f = open(new_filename, 'w')

    for time_used in time_used_list:
        f.write(str(time_used))
        f.write("\n")


def write_cuts(header, filename, nb_iter, nb_cuts, time_used_list, 
               expanded_list, generated_list, cuts_per_iter, sp_time_used_list, mp_cost_list):
    sep = filename.split('/')
    print(sep)

    folder = "stats/" + header + "/" + sep[2] + "/" + sep[3] + "/"
    os.makedirs(folder, exist_ok=True)

    new_filename = "stats/" + header + "/" + sep[2] + "/" + sep[3] + "/" + sep[4][:-3] + "times"
    print(new_filename)

    f = open(new_filename, 'w')

    f.write(str(nb_iter) + ' ')
    f.write(str(nb_cuts) + '\n')

    for i in range(len(time_used_list)):
        f.write(str(time_used_list[i]) + ' ')
        f.write(str(sp_time_used_list[i]) + ' ')
        f.write(str(expanded_list[i]) + ' ')
        f.write(str(generated_list[i]) + ' ')
        f.write(str(cuts_per_iter[i]) + ' ')
        f.write(str(mp_cost_list[i]) + ' ')
        f.write("\n")


def read_xls():
    filename = "albp2.xlsx"
    df = pd.read_excel(filename)
    names = df['name'].to_list()
    cs = df['c'].to_list()
    ms = df['m'].to_list()

    mlist = {}
    for i in range(len(names)):
        s = names[i]
        mlist[s.lower() + '_c=' + str(cs[i]) + '.alb'] = ms[i]

    return mlist



def get_ub():
    filename = "albp2.xlsx"
    df = pd.read_excel(filename)
    names = df['name'].to_list()
    ubs = df['ub'].to_list()
    cs = df['c'].to_list()
    alphas = df['alpha'].to_list()

    ublist = {}
    for i in range(len(names)):
        s = names[i]
        alpha = alphas[i]
        if alpha == 0.5:
            ublist["0.50" + '/' + s.lower() + '_c=' + str(cs[i]) + '.alb'] = ubs[i]
        elif alpha == 1.0:
            ublist["1.00" + '/' + s.lower() + '_c=' + str(cs[i]) + '.alb'] = ubs[i]
        else:
            ublist[str(alpha) + '/' + s.lower() + '_c=' + str(cs[i]) + '.alb'] = ubs[i]

    return ublist




def get_lb():
    filename = "albp2.xlsx"
    df = pd.read_excel(filename)
    names = df['name'].to_list()
    lbs = df['lb'].to_list()
    cs = df['c'].to_list()
    alphas = df['alpha'].to_list()

    lblist = {}
    for i in range(len(names)):
        s = names[i]
        alpha = alphas[i]
        if alpha == 0.5:
            lblist["0.50" + '/' + s.lower() + '_c=' + str(cs[i]) + '.alb'] = lbs[i]
        elif alpha == 1.0:
            lblist["1.00" + '/' + s.lower() + '_c=' + str(cs[i]) + '.alb'] = lbs[i]
        else:
            lblist[str(alpha) + '/' + s.lower() + '_c=' + str(cs[i]) + '.alb'] = lbs[i]

    return lblist




def read_obj_at_time(filename, m):
    print("filename in read_obj_at_time:", filename)
    obj_at_time = []
    with open(filename) as f:
        line = f.readline().strip().split()
        if line[0] == 'None':
            return obj_at_time
        for i in range(m):
            line = f.readline()
        line = f.readline().strip()
        if len(line.split()) == 0:
            line = f.readline().strip()
        while line:
            params = line.split()
            if len(line) > 0:
                obj_at_time.append((float(params[0]), int(float(params[1]))))
            line = f.readline().strip()
    return obj_at_time



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


# merge the above two functions
def compute_earliest_latest(number_of_tasks, number_of_stations, 
                            ub_cycle, task_times, predecessors, followers,
                            forward_min_task_f, backward_min_task_f,
                            forward_min_task_b, backward_min_task_b,
):
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



def calculate_Ai(number_of_tasks, earliest, latest):
    Ai = {}
    for i in range(number_of_tasks):
        Ai[i+1] = []
        for j in range(number_of_tasks):
            if latest[i+1] < earliest[j+1] or latest[j+1] < earliest[i+1]:
                Ai[i+1].append(j+1)
    return Ai





def calculate_F_and_P(number_of_tasks, task_times, 
    predecessors, followers,
    all_predecessors, all_followers, 
    forward, backward,
    Ai):

    FiF = {}
    PiF = {}

    F_diff = {}
    for i in range(number_of_tasks):
        F_diff[i+1] = []
        for j in range(number_of_tasks):
            if j+1 in all_followers[i+1] and j+1 not in followers[i+1]:
                F_diff[i+1].append(j+1)

    for i in range(number_of_tasks):
        FiF[i+1] = []
        for j in range(number_of_tasks):
            if j+1 not in F_diff[i+1] and \
                j+1 not in all_predecessors[i+1] and \
                j+1 not in Ai[i+1] and \
                j+1 != i+1:
                FiF[i+1].append(j+1)
    
    for i in range(number_of_tasks):
        PiF[i+1] = []
        for j in range(number_of_tasks):
            if i+1 in FiF[i+1]:
                PiF[i+1].append(j+1)

    FiB = {}
    PiB = {}

    for i in range(number_of_tasks):
        FiB[i+1] = []
        for j in range(number_of_tasks):
            if j+1 not in all_followers[i+1] and \
                j+1 not in Ai[i+1]:
                FiB[i+1].append(j+1)
    
    for i in range(number_of_tasks):
        PiB[i+1] = []
        for j in range(number_of_tasks):
            if i+1 in FiB[i+1]:
                PiB[i+1].append(j+1)

    return FiF, PiF, FiB, PiB





def calculate_pair_earliest(number_of_tasks, task_times, cycle_time, m_ub,
    all_predecessors, all_followers, forward, backward,
    FiF, PiF, FiB, PiB):

    aijF = {}
    for i in range(number_of_tasks):
        for j in FiF[i+1]:
            aijF_ = task_times[i+1] + task_times[j] + forward[i][j-1]
            for k in range(number_of_tasks):
                if (k+1 in all_predecessors[i+1] or k+1 in all_predecessors[j]) and \
                    k+1!=i+1 and k+1!=j:
                    aijF_ += task_times[k+1]
            aijF_ = math.ceil(aijF_ / cycle_time)
            aijF[(i+1,j)] = aijF_
    
    bijF = {}
    for i in range(number_of_tasks):
        for j in FiF[i+1]:
            bijF_ = task_times[i+1] + task_times[j] + forward[i][j-1]
            for k in range(number_of_tasks):
                if (k+1 in all_followers[i+1] or k+1 in all_followers[j]) and \
                    k+1!=i+1 and k+1!=j:
                    bijF_ += task_times[k+1]
            bijF_ = m_ub + 1 - math.ceil(bijF_ / cycle_time)
            bijF[(i+1,j)] = bijF_

    aijB = {}
    for i in range(number_of_tasks):
        for j in FiB[i+1]:
            aijB_ = task_times[i+1] + task_times[j] + backward[i][j-1]
            for k in range(number_of_tasks):
                if (k+1 in all_predecessors[i+1] or k+1 in all_predecessors[j]) and \
                    k+1!=i+1 and k+1!=j:
                    aijB_ += task_times[k+1]
            aijB_ = math.ceil(aijB_ / cycle_time)
            aijB[(i+1,j)] = aijB_

    bijB = {}
    for i in range(number_of_tasks):
        for j in FiB[i+1]:
            bijB_ = task_times[i+1] + task_times[j] + backward[i][j-1]
            for k in range(number_of_tasks):
                if (k+1 in all_followers[i+1] or k+1 in all_followers[j]) and \
                    k+1!=i+1 and k+1!=j:
                    bijB_ += task_times[k+1]
            bijB_ = m_ub + 1 - math.ceil(bijB_ / cycle_time)
            bijB[(i+1,j)] = bijB_

    return aijF, bijF, aijB, bijB





def calculate_B_and_R(number_of_tasks, earliest, latest,
                      all_predecessors, all_followers):
    Bi = {}
    for i in range(number_of_tasks):
        Bi[i+1] = []
        for j in range(number_of_tasks):
            if latest[i+1] < earliest[j+1] or latest[j+1] < earliest[i+1]:
                Bi[i+1].append(j+1)

    Rij = {}
    for i in range(number_of_tasks):
        for j in all_followers[i+1]:
            if j not in Bi[i+1]:
                Rij[(i+1,j)] = []
                for k in range(number_of_tasks):
                    if k+1 in all_followers[i+1] and k+1 in all_predecessors[j]:
                        Rij[(i+1,j)].append(k+1)
                if i+1 not in Rij[(i+1,j)]:
                    Rij[(i+1,j)].append(i+1)
                if j not in Rij[(i+1,j)]:
                    Rij[(i+1,j)].append(j)

    return Bi, Rij




def calculate_delta(number_of_tasks, earliest, latest,
                    all_predecessors, all_followers,
                    forward, backward,
                    Bi, Rij,
                    FiF, PiF, FiB, PiB):
    delta = {}
    for i in range(number_of_tasks):
        for j in all_followers[i+1]:
            if j not in Bi[i+1]:
                for v in Rij[(i+1, j)]:
                    if v != j:
                        # too long for loop
                        delta_ = 1000000
                        for vp in Rij[(i+1, j)]:
                            if vp in FiF[v]:
                                delta_ = min(delta_, forward[v-1][vp-1])
                        delta[(i+1,v,j)] = delta_

    deltap = {}
    for i in range(number_of_tasks):
        for j in all_followers[i+1]:
            if j not in Bi[i+1]:
                for v in Rij[(i+1, j)]:
                    if v != i+1:
                        # too long for loop
                        delta_ = 1000000
                        for vp in Rij[(i+1, j)]:
                            if vp in PiF[v]:
                                delta_ = min(delta_, forward[vp-1][v-1])
                        deltap[(i+1,v,j)] = delta_
    
    return delta, deltap



def calculate_dij(number_of_tasks, task_times, cycle_time,
                earliest, latest,
                all_predecessors, all_followers,
                forward, backward,
                Bi, Rij,
                FiF, PiF, FiB, PiB,
                delta, deltap):
    dij = {}
    for i in range(number_of_tasks):
        for j in all_followers[i+1]:
            if j not in Bi[i+1]:
                dij_ = backward[j-1][i]
                for v in Rij[(i+1, j)]:
                    dij_ += task_times[v]
                t1 = 0
                for v in Rij[(i+1, j)]:
                    if v != j:
                        t1 += delta[(i+1, v, j)]
                t2 = 0
                for v in Rij[(i+1, j)]:
                    if v != i+1:
                        t2 += deltap[(i+1, v, j)]
                dij_ += max(t1, t2)
                dij[(i+1,j)] = dij_
                dij[(j,i+1)] = dij_

    Ai = Bi.copy()
    for i in range(number_of_tasks):
        for j in range(number_of_tasks):
            if (i+1, j+1) in dij and dij[(i+1, j+1)] > cycle_time and j+1 not in Ai[i+1]:
                Ai[i+1].append(j+1)
    
    return dij, Ai




def calculate_gamma(number_of_tasks, forward, backward,
                    FiF, PiF, FiB, PiB):
    gamma_i1 = {}
    for i in range(number_of_tasks):
        option1 = 1000000
        for j in FiF[i+1]:
            option1 = min(option1, forward[i][j-1])
        option2 = 1000000
        for j in FiB[i+1]:
            option2 = min(option2, backward[i][j-1])
        gamma_i1[i+1] = min(option1, option2)
    
    gamma_i2 = {}
    for i in range(number_of_tasks):
        option1 = 1000000
        for j in PiF[i+1]:
            option1 = min(option1, forward[j-1][i]-gamma_i1[i+1])
        option2 = 1000000
        for j in PiB[i+1]:
            option2 = min(option2, backward[j-1][i]-gamma_i1[i+1])
        gamma_i2[i+1] = min(option1, option2)

    sdF = []
    for i in range(number_of_tasks):
        for j in range(number_of_tasks):
            sdF.append(forward[i][j] - gamma_i1[i+1] - gamma_i2[j+1])
    sdF_sorted = sorted(sdF)

    sdB = []
    for i in range(number_of_tasks):
        for j in range(number_of_tasks):
            sdB.append(backward[i][j] - gamma_i1[i+1] - gamma_i2[j+1])
    sdB_sorted = sorted(sdB)

    return gamma_i1, gamma_i2, sdF_sorted, sdB_sorted




def calculate_lb_cycle_new(number_of_tasks, number_of_stations, task_times,
    forward, backward, gamma_i1, gamma_i2, sdF_sorted, sdB_sorted):
    option1 = 0
    for i in range(number_of_tasks):
        option1 = max(option1, task_times[i+1]+backward[i][i])
    option2 = 0
    for i in range(number_of_tasks):
        option2 += (task_times[i+1] + gamma_i1[i+1] + gamma_i2[i+1])
    for s in range(number_of_tasks-number_of_stations):
        option2 += sdF_sorted[s]
    for s in range(number_of_stations):
        option2 += sdB_sorted[s]
    option2 = math.ceil(option2 / number_of_stations)

    return max(option1, option2)



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



def compute_ub_cycle_time_new(
    number_of_tasks, number_of_stations, task_times,
    all_predecessors, all_followers,
    forward, backward, gamma_i1, gamma_i2,
    FiF, PiF, FiB, PiB
):
    '''
    Compute an upper bound of cycle time
    '''
    tau_p = []
    for i in range(number_of_tasks):
        tau_p_ = 0
        for j in FiF[i+1]:
            tau_p_ = max(tau_p_, forward[i][j-1])
        tau_p.append(tau_p_)
    tau_pp = []
    for i in range(number_of_tasks):
        tau_pp_ = 0
        for j in PiF[i+1]:
            tau_pp_ = max(tau_pp_, forward[j-1][i])
        tau_pp.append(tau_pp_)
    mu_p = []
    for i in range(number_of_tasks):
        mu_p_ = 0
        for j in FiB[i+1]:
            mu_p_ = max(mu_p_, backward[i][j-1])
        mu_p.append(mu_p_)
    mu_pp = []
    for i in range(number_of_tasks):
        mu_pp_ = 0
        for j in PiB[i+1]:
            mu_pp_ = max(mu_pp_, backward[j-1][i])
        mu_pp.append(mu_pp_)
    
    tau_p_sorted = sorted(tau_p, reverse=True)
    tau_pp_sorted = sorted(tau_pp, reverse=True)
    mu_p_sorted = sorted(mu_p, reverse=True)
    mu_pp_sorted = sorted(mu_pp, reverse=True)

    t1 = 0
    t2 = 0
    for i in range(number_of_tasks-number_of_stations):
        t1 += tau_p_sorted[i]
        t2 += tau_pp_sorted[i]
    tau_hat = min(t1, t2)

    m1 = 0
    m2 = 0
    for i in range(number_of_stations):
        m1 += mu_p_sorted[i]
        m2 += mu_pp_sorted[i]
    mu_hat = min(m1, m2)

    
    option1 = 0
    for i in range(number_of_tasks):
        if task_times[i+1] + forward[i][i] > option1:
            option1 = task_times[i+1] + forward[i][i]

    option2 = 0
    for i in range(number_of_tasks):
        option2 += task_times[i+1]
    
    if number_of_stations // 2 == 0:
        option2 = math.floor(
            2*(option2 + tau_hat + mu_hat -1) / number_of_stations
        )
    else:
        option2 = math.floor(
            2*(option2 + tau_hat + mu_hat) / (number_of_stations+1)
        )

    return max(option1, option2)




def calculate_new_earliest_and_latest(
    number_of_tasks, number_of_stations, task_times, ub_cycle,
    all_predecessors, all_followers,
    forward, backward, gamma_i1, gamma_i2):

    earliestn = {}
    for i in range(number_of_tasks):
        e = task_times[i+1] + gamma_i1[i+1] + gamma_i2[i+1]
        for j in all_predecessors[i+1]:
            e += (task_times[j] + gamma_i1[j] + gamma_i2[j])
        e = math.ceil(e / ub_cycle)
        earliestn[i+1] = e
    
    latestn = {}
    for i in range(number_of_tasks):
        e = task_times[i+1] + gamma_i1[i+1] + gamma_i2[i+1]
        for j in all_followers[i+1]:
            e += (task_times[j] + gamma_i1[j] + gamma_i2[j])
        e = number_of_stations + 1 - math.ceil(e / ub_cycle)
        latestn[i+1] = e

    return earliestn, latestn


# lower and upper bound of the number of stations
def compute_m_bounds(number_of_tasks, cycle_time, task_times,
                     forward_min_task, backward_min_task):
    lb = math.ceil(
        (sum(task_times.values()) + 
        sum([min(forward_min_task[i], backward_min_task[i]) for i in range(number_of_tasks)])
        ) / cycle_time
    )
    ub = min(2 * lb, number_of_tasks)
    return lb, ub


# compute the lower bound of the cycle time
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


# compute the contribution of a task to the station time
def compute_contribution(
    number_of_tasks, number_of_stations, task_times,
    forward_max_task_f, backward_max_task_f,
    forward_max_task_b, backward_max_task_b,
):
    contribution = {}
    for i in range(number_of_tasks):
        contri = 0
        contri += task_times[i+1]
        t1 = backward_max_task_b[i] + forward_max_task_f[i]
        t2 = forward_max_task_b[i] + forward_max_task_f[i]
        t3 = forward_max_task_b[i] + backward_max_task_f[i]
        contri += max(t1, t2, t3)
        contribution[i+1] = contri

    return contribution


# ! compute the contribution corresponding to the specific subproblem
def compute_contribution_sp(
    tasks_in_sp, task_times, contribution,
    FiF, PiF, FiB, PiB, forward, backward,  
):
    contribution_item = {}
    for i in tasks_in_sp:
        contri = 0
        contri += task_times[i]
        # t1
        t1_left = 0
        for j in tasks_in_sp:
            if j in FiF[i] and j != i:
                t1_left = max(t1_left, forward[i-1][j-1])
        t1_right = 0
        for j in tasks_in_sp:
            if j in PiB[i] and j != i:
                t1_right = max(t1_right, backward[j-1][i-1])
        # t2
        t2_left = 0
        for j in tasks_in_sp:
            if j in PiF[i] and j != i:
                t2_left = max(t2_left, forward[j-1][i-1])
        t2_right = 0
        for j in tasks_in_sp:
            if j in FiB[i] and j != i:
                t2_right = max(t2_right, backward[i-1][j-1])
        # t3
        t3_left = 0
        for j in tasks_in_sp:
            if j in PiF[i] and j != i:
                t3_left = max(t3_left, forward[j-1][i-1])
        t3_right = 0
        for j in tasks_in_sp:
            if j in FiF[i] and j != i:
                t3_right = max(t3_right, forward[i-1][j-1])
        contri += max(t1_left + t1_right, 
                      t2_left + t2_right, 
                      t3_left + t3_right)
        contribution_item[i] = contri
    contribution.append(contribution_item)


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