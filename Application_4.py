"""comparison of proteins"""
import random
from matplotlib import pyplot as plt


def load_protein(file_name):
    with open(file_name, 'r') as f:
        return f.readline()


def load_scoring_matrix(file_name):
    lines = []
    with open(file_name, 'r') as f:
        for line in f.readlines():
            lines.append(line.split())

    matrix = {}
    for letter in lines[0]:
        matrix[letter] = {}

    first_line = lines[0]
    lines = lines[1:]

    for line in lines:
        for idx in range(len(first_line)):
            matrix[line[0]][first_line[idx]] = int(line[idx+1])

    return matrix

human_eyeless = load_protein('alg_HumanEyelessProtein.txt')
fruitfly_eyeless = load_protein('alg_FruitflyEyelessProtein.txt')
PAM50 = load_scoring_matrix('alg_PAM50.txt')


def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
    """
    Creates scoring matrix based on input.
    :param alphabet: set of letters
    :param diag_score: score when letters are same
    :param off_diag_score: score when letters are different
    :param dash_score: score when letter matched with dash
    :return: dictionary of dictionaries
    """
    matrix = {'-': {'-': dash_score}}

    for letter in alphabet:
        interior = {}
        for other in alphabet:
            if other == letter:
                interior[other] = diag_score
            else:
                interior[other] = off_diag_score
        interior['-'] = dash_score
        matrix[letter] = interior
        matrix['-'][letter] = dash_score

    return matrix


def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag):
    """
    Computes alignment scores for input sequences.
    :param seq_x: first sequence
    :param seq_y: second sequence
    :param scoring_matrix: scoring matrix
    :param global_flag: if true then global alignment, if false then local alignment
    :return: 2d matrix grid of maximum score for each alignment
    """
    x_len = len(seq_x)
    y_len = len(seq_y)
    matrix = [[None for dummyidx in range(y_len+1)] for dummyidx in range(x_len+1)]

    if len(matrix) == 0:
        print "Empty"
        return [[0]]

    matrix[0][0] = 0

    for pos_i in range(1, x_len+1):
        temp_score = 0 if matrix[pos_i-1][0] + scoring_matrix[seq_x[pos_i-1]]['-'] < 0 and not global_flag \
            else matrix[pos_i-1][0] + scoring_matrix[seq_x[pos_i-1]]['-']
        matrix[pos_i][0] = temp_score

    for pos_j in range(1, y_len+1):
        temp_score = 0 if matrix[0][pos_j-1] + scoring_matrix['-'][seq_y[pos_j-1]] < 0 and not global_flag \
            else matrix[0][pos_j-1] + scoring_matrix['-'][seq_y[pos_j-1]]
        matrix[0][pos_j] = temp_score

    for pos_i in range(1, x_len+1):
        for pos_j in range(1, y_len+1):
            opt1 = matrix[pos_i-1][pos_j-1] + scoring_matrix[seq_x[pos_i-1]][seq_y[pos_j-1]]
            opt2 = matrix[pos_i-1][pos_j] + scoring_matrix[seq_x[pos_i-1]]['-']
            opt3 = matrix[pos_i][pos_j-1] + scoring_matrix['-'][seq_y[pos_j-1]]
            temp_score = 0 if max(opt1, opt2, opt3) < 0 and not global_flag else max(opt1, opt2, opt3)
            matrix[pos_i][pos_j] = temp_score

    return matrix


def compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    Computes optimal global alignment based on alignment and scoring matrices
    :param seq_x: first sequence
    :param seq_y: second sequence
    :param scoring_matrix: scoring matrix
    :param alignment_matrix: alignment matrix
    :return: tuple (score, x sequence, y sequence)
    """
    x_len, y_len = len(seq_x), len(seq_y)
    x_prime, y_prime = '', ''

    while x_len != 0 and y_len != 0:
        if alignment_matrix[x_len][y_len] == alignment_matrix[x_len-1][y_len-1] + \
                scoring_matrix[seq_x[x_len-1]][seq_y[y_len-1]]:
            x_prime = seq_x[x_len-1] + x_prime
            y_prime = seq_y[y_len-1] + y_prime
            x_len -= 1
            y_len -= 1
        else:
            if alignment_matrix[x_len][y_len] == alignment_matrix[x_len-1][y_len] + scoring_matrix[seq_x[x_len-1]]['-']:
                x_prime = seq_x[x_len-1] + x_prime
                y_prime = '-' + y_prime
                x_len -= 1
            else:
                x_prime = '-' + x_prime
                y_prime = seq_y[y_len-1] + y_prime
                y_len -= 1

    while x_len != 0:
        x_prime = seq_x[x_len-1] + x_prime
        y_prime = '-' + y_prime
        x_len -= 1

    while y_len != 0:
        x_prime = '-' + x_prime
        y_prime = seq_y[y_len-1] + y_prime
        y_len -= 1

    score = 0
    for idx in range(len(x_prime)):
        score += scoring_matrix[x_prime[idx]][y_prime[idx]]

    return tuple([score, x_prime, y_prime])


def compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    computes optimale local alignment based on scoring and alignment matrices
    :param seq_x: first sequence
    :param seq_y: second sequence
    :param scoring_matrix: scoring matrix
    :param alignment_matrix: alignment matrix
    :return: tuple (score, x sequence, y sequence)
    """
    x_len, y_len = len(seq_x), len(seq_y)
    highest = alignment_matrix[x_len][y_len]
    for idx in range(x_len+1)[::-1]:
        for idx2 in range(y_len+1)[::-1]:
            if alignment_matrix[idx][idx2] > highest:
                highest = alignment_matrix[idx][idx2]
                x_len = idx
                y_len = idx2

    x_prime, y_prime = '', ''

    while x_len != 0 and y_len != 0:
        if alignment_matrix[x_len][y_len] == 0:
            break

        if alignment_matrix[x_len][y_len] == alignment_matrix[x_len-1][y_len-1] + \
                scoring_matrix[seq_x[x_len-1]][seq_y[y_len-1]]:
            x_prime = seq_x[x_len-1] + x_prime
            y_prime = seq_y[y_len-1] + y_prime
            x_len -= 1
            y_len -= 1
        else:
            if alignment_matrix[x_len][y_len] == alignment_matrix[x_len-1][y_len] + scoring_matrix[seq_x[x_len-1]]['-']:
                x_prime = seq_x[x_len-1] + x_prime
                y_prime = '-' + y_prime
                x_len -= 1
            else:
                x_prime = '-' + x_prime
                y_prime = seq_y[y_len-1] + y_prime
                y_len -= 1

    score = 0
    for idx in range(len(x_prime)):
        score += scoring_matrix[x_prime[idx]][y_prime[idx]]

    return tuple([score, x_prime, y_prime])

# alignment = compute_alignment_matrix(human_eyeless, fruitfly_eyeless, PAM50, False)
# local_alignments = compute_local_alignment(human_eyeless, fruitfly_eyeless, PAM50, alignment)


consensus_PAX_domain = load_protein('alg_ConsensusPAXDomain.txt')
# human_alignment = local_alignments[1].replace('-', '')
# h_alignment_matrix = compute_alignment_matrix(human_alignment, consensus_PAX_domain, PAM50, True)
# fruitfly_alignment = local_alignments[2].replace('-', '')
# f_alignment_matrix = compute_alignment_matrix(fruitfly_alignment, consensus_PAX_domain, PAM50, True)
#
# h_consensus_alignment = compute_global_alignment(human_alignment, consensus_PAX_domain, PAM50, h_alignment_matrix)
# f_consensus_alignment = compute_global_alignment(fruitfly_alignment, consensus_PAX_domain, PAM50, f_alignment_matrix)


def compare_similarity(str1, str2):
    total = len(str1)
    similar = 0.0
    for idx in range(total):
        if str1[idx] == str2[idx]:
            similar += 1

    return similar / total * 100

# print len(f_consensus_alignment[1])
# print compare_similarity(h_consensus_alignment[1], h_consensus_alignment[2]) * len(h_consensus_alignment[1]) / 100
# print compare_similarity(f_consensus_alignment[1], f_consensus_alignment[2]) * len(f_consensus_alignment[1]) / 100


def shuffle_string(string):
    split = list(string)
    random.shuffle(split)
    return ''.join(split)


def generate_null_distribution(seq_x, seq_y, scoring_matrix, num_trials):
    scoring_distribution = {}
    for trial in range(num_trials):
        print "Doing trial number ", trial + 1
        rand_y = shuffle_string(seq_y)
        score = compute_local_alignment(seq_x, rand_y, scoring_matrix,
                                        compute_alignment_matrix(seq_x,rand_y, scoring_matrix, False))[0]
        if score in scoring_distribution:
            scoring_distribution[score] += 1
        else:
            scoring_distribution[score] = 1

    return scoring_distribution

# null_human_fruitfly = generate_null_distribution(human_eyeless, fruitfly_eyeless, PAM50, 1000)
# print "Done generating distribution"
# print null_human_fruitfly
#
# normalized_n_h_f = {key: null_human_fruitfly[key] / 1000.0 for key in null_human_fruitfly}
# print normalized_n_h_f


def list_from_norm_dict(dict):
    full = []
    for key in dict:
        full.extend([key]*int(dict[key]*1000))

    return full


def mean(data):
    length = float(len(data))
    return sum(data)/length


def _ss(data):
    avg = mean(data)
    squares = sum((point - avg)**2 for point in data)
    return squares


def stdev(data, ddof=0):
    length = float(len(data))
    squares = _ss(data)
    pvar = squares/(length-ddof)
    return pvar**0.5


def z_score(score, average, standard_deviation):
    return (score - average) / standard_deviation


# dist_list = list_from_norm_dict(normalized_n_h_f)
# print dist_list
#
# avg = mean(dist_list)
# sigma = stdev(dist_list)
#
# print "Mean is ", avg
# print "Standard deviation is ", sigma
# print "Z-score is ", z_score(compute_local_alignment(human_eyeless, fruitfly_eyeless, PAM50,
#                                                      compute_alignment_matrix(human_eyeless, fruitfly_eyeless,
#                                                                               PAM50, False))[0],
#                              avg, sigma)
#
#
# plt.bar([key for key in normalized_n_h_f], [normalized_n_h_f[key] for key in normalized_n_h_f])
# plt.ylabel('Proportion of trials')
# plt.xlabel('Score')
# plt.title('Normalized Null Distribution of Human and Fruitfly Eyeless Proteins')
# plt.show()


def edit_distance(seq_x, seq_y, scoring_matrix, global_flag=True):
    score = compute_global_alignment(seq_x, seq_y, scoring_matrix,
                                     compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag))

    return len(seq_x) + len(seq_y) - score[0]


def load_word_list(file_name):
    words = []
    with open(file_name, 'r') as f:
        for line in f.readlines():
            words.append(line.rstrip())

    return words

dictionary = load_word_list('assets_scrabble_words3.txt')


def check_spelling(checked_word, dist, word_list):
    scoring = build_scoring_matrix({'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q',
                                    'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'}, 2, 1, 0)
    possible = set([])
    for word in word_list:
        if edit_distance(checked_word, word, scoring) <= dist:
            possible.add(word)

    return possible

print check_spelling('humble', 1, dictionary)
print check_spelling('firefly', 2, dictionary)

