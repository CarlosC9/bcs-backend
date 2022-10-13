def get_marginal_l_estimate(filename: str):
    with open(filename, 'r') as f:
        for l in f:
            if l.strip().startswith('marginal L estimate ='):
                return float(l.split("=")[1])
    return "error"

if __name__ == '__main__':
    import sys
    birthDeathModel = sys.argv[1]
    yuleModel = sys.argv[2]

    ml_bd = get_marginal_l_estimate(birthDeathModel)
    ml_y = get_marginal_l_estimate(yuleModel)

    better = "Birth Death Model"
    worse = "Yule Model"
    better_n = ml_bd
    worse_n = ml_y
    if ml_y > ml_bd:
        better = "Yule Model"
        worse = "Birth Death Model"
        better_n = ml_y
        worse_n = ml_bd
    bf = abs(worse_n) - abs(better_n)#it must be in the opposite way because they are negative numbers
    if bf == 0:
        output_str = f"Bayes Factor = {bf}.\n{better} and {worse} are equally valid."
    elif bf <= 1.1: 
        output_str = f"Bayes Factor = {bf}.\nThe {better} gets slight support over {worse}."
    elif bf <= 3:
        output_str = f"Bayes Factor = {bf}.\n{better} gets positive support over {worse}."
    elif bf <= 5:
        output_str = f"Bayes Factor = {bf}.\n{better} gets strong positive support over {worse}."
    else:
        output_str = f"Bayes Factor = {bf}.\n{better} gets overwhelming positive support over {worse}."

    with open("bayes_factor.txt", "w") as f:
        f.write(output_str)
