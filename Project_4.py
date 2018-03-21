"""implementation of project 4 for algorithmic thinking"""


def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
    """
    Creates scoring matrix based on input.
    :param alphabet: set of letters
    :param diag_score: score when letters are same
    :param off_diag_score: score when letters are different
    :param dash_score: score when letter matched with dash
    :return: dictionary of dictionaries
    """
    matrix = {}
    matrix['-'] = {'-': dash_score}

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
        if alignment_matrix[x_len][y_len] == alignment_matrix[x_len-1][y_len-1] + scoring_matrix[seq_x[x_len-1]][seq_y[y_len-1]]:
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

    while  y_len != 0:
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

    print highest, x_len, y_len, alignment_matrix
    x_prime, y_prime = '', ''

    while x_len != 0 and y_len != 0:
        if alignment_matrix[x_len][y_len] == 0:
            break

        if alignment_matrix[x_len][y_len] == alignment_matrix[x_len-1][y_len-1] + scoring_matrix[seq_x[x_len-1]][seq_y[y_len-1]]:
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

