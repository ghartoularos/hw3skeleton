# coding: utf-8

# Smith-Waterman Algorithm
'''

Implement the Smith-Waterman algorithm and instrument the code such that it can 
use any scoring matrix provided (i.e. will read it in from a separate file). You 
may adapt code you can obtain from the Web (there are a number implementations 
you can find). But you must demonstrate your understanding of the algorithm with 
detailed comments.

The Smith-Waterman algorithm is an algorithm for determining similar regions 
between two strings of nucleic acid sequences or protein sequences. Instead of 
looking at the entire sequence, the Smith–Waterman algorithm compares segments of 
all possible lengths and optimizes the similarity measure.

'''
# Smith-Waterman Algorithm: https://gist.github.com/radaniba/11019717

# Import packages
import pandas as pd


def get_submat(submatfile):
    '''
    Reads in a substitution matrix file as a pandas dataframe
    '''
    with open(submatfile,'r') as f:
        count = sum([i[0] == '#' for i in f.read().strip().split('\n')])
    matrix = pd.read_csv(submatfile,sep=' ',skiprows=count,skipinitialspace=True,index_col=False)
    matrix.index = matrix.columns
    return matrix

def main(submatfile,seq1,seq2,gapopen=-11,gapext=-1,road='high',q1=False):
    '''
    Main function calls all other functions and prints results
    '''

    assert gapopen <= 0 and gapext <= 0, 'Gap penalties must be non-positive.'
    assert type(seq1) == str and type(seq2) == str, 'Sequences are not strings.'
    # Bring in the substitution matrix file as submat
    try:
        submat = get_submat(submatfile)
    except:
        print('Could not retrieve substitution matrix.')
        return
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    # The scoring matrix contains an extra row and column to allow for 
    # end gap alignments.
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    # Create the scoring matrix
    scoremat, dirmat, start_pos = create_scoremat(submat, rows, cols, seq1, 
                                                  seq2, gapopen, gapext, road, q1)
    if q1:
        max_score = scoremat # Maxscore was passed as scoremat
        return max_score
    # Traceback. Find the optimal path through the scoring matrix. This path
    # corresponds to the optimal local sequence alignment.
    seq1_aligned, seq2_aligned = traceback(seq1, seq2, scoremat, dirmat, start_pos)
    
    # The aligned sequences should be the same size, since if a position from either
    # sequence was gapped, it was replaced with a '-' character.
    assert len(seq1_aligned) == len(seq2_aligned), 'Aligned sequences are not the same size.'

    # Pretty print the results. The printing follows the format of BLAST results
    # as closely as possible.
    alignment_str, idents, gaps, mismatches = alignment_string(seq1_aligned, seq2_aligned)
    alength = len(seq1_aligned)
    print('')
    print(' Identities = {0}/{1} ({2:.1%}), Gaps = {3}/{4} ({5:.1%})'.format(idents,
          alength, idents / float(alength), gaps, alength, gaps / float(alength)))
    print('')
    for i in range(0, alength, 60):
        seq1_slice = seq1_aligned[i:i+60]
        print('Seq1:  {0:<4}  {1}  {2:<4}'.format(i + 1, seq1_slice, i + len(seq1_slice)))
        print('             {0}'.format(alignment_str[i:i+60]))
        seq2_slice = seq2_aligned[i:i+60]
        print('Seq2:  {0:<4}  {1}  {2:<4}'.format(i + 1, seq2_slice, i + len(seq2_slice)))
        print('')

def create_scoremat(submat, rows, cols, seq1, seq2, 
                    gapopen, gapext, road, q1):
    '''
    Create a matrix of scores representing trial alignments of the two sequences.
    Sequence alignment can be treated as a graph search problem. This function
    creates a graph (2D matrix) of scores, which are based on trial alignments
    of different base pairs. The path with the highest cummulative score is the
    best alignment.
    
    Some things about this implementation:
    
    (1) It uses an affine gap penalty scheme: the penalty for opening a gap is some constant
    while the penalty for extending a gap is linear with the length of the gap 
    Other implementations can use a simple constant, linear, or convex gap penalty scheme, 
    which all perform differently.
    
    (2) It uses the first max score encountered as the start position of the traceback.
    This algorithm can theoretically find multiple optimal alignements, but since
    they will all have the same score, deciding which one is better is beyond the scope
    of the algorithm. 
    
    (3) It has the option for "road" to choose the "high road" or "low road" when there
    are multiple alignments yielded the max score. High road check for max score in the 
    high -> mid -> low order, whereas low road checks in the low -> mid -> high order.
    '''
    # Scoring matrix is initialized to zeros. Many zeros in the inner matrix
    # (excluding first row, first column) will be replaced with score. First
    # row and column will remain zero for end gaps, which are not penalized.
    scoremat = [[0 for col in range(cols)] for row in range(rows)]
    
    gapup, gapleft = gapopen, gapopen
    # Keep track of where the max score came from in a "direction matrix"
    # so that traceback can easily be done.
    dirmat = [[0 for col in range(cols)] for row in range(rows)]
    # Fill the scoring matrix.
    max_score = 0
    max_pos   = None    # The row and column of the highest score in matrix.
    for i in range(1, rows):
        for j in range(1, cols):
            # Get the score for a certain entry in the matrix
            score, dirmat = calc_score(seq1, seq2, scoremat, submat, dirmat, 
                                       gapup, gapleft, i, j, road)
            if dirmat[i][j] == 3:
                gapup = gapext
            elif dirmat[i][j] == 1:
                gapleft = gapext
            else:
                gapup, gapleft = gapopen, gapopen
            # The max score is where we will start the traceback
            if score > max_score:
                max_score = score
                max_pos   = (i, j)
            # Fill in the scoremat with the score found for the given position
            scoremat[i][j] = score
    # In case, for some reason, you didn't find a score > 0:
    assert max_pos is not None, 'The x, y position with the highest score was not found.'
    if q1:
        return max_score, None, None
    return scoremat, dirmat, max_pos

def calc_score(seq1, seq2, scoremat, submat, dirmat, gapup, 
               gapleft, x, y, road='high'):
    '''
    Calculate score for a given x, y position in the scoring matrix.
    The score is based on the up, left, and upper-left neighbors.
    '''
    # The amino acid in question from each sequence is the row or column position of the
    # score matrix *minus 1* since the scoremat has one more zero-initialized
    # row and column at the start to account for end gaps
    n1 = seq1[x - 1] 
    n2 = seq2[y - 1]
    
    # Similarity score from the substitution matrix
    similarity = submat[n1][n2]
    
    # The score derived from the diagonal will take the amino acid 
    # from both sequences, regardless of whether or not they match
    diag_score = scoremat[x - 1][y - 1] + similarity
    
    # The up and left scores take an amino acid from only 
    # one sequence and include a gap in the other
    up_score   = scoremat[x - 1][y] + gapup
    left_score = scoremat[x][y - 1] + gapleft
    
    # Returning the max of all the scores *and zero* is what makes this a 
    # local alignment as opposed to a global alignment
    scores = (0, left_score, diag_score, up_score)
    
    # Keep track of where the score came from
    maxscore, dirmat = find_move(scores, dirmat, road, x, y)

    return maxscore, dirmat

def find_move(scores, dirmat, road, x, y):
    '''
    Source of the score: symbol
    
    Up: 3
    Diagonal: 2
    Left: 1
    No high score, end gap: 0
    
    In the case that the score could have come from multiple adjacent cells,
    there are multiple possible alignments. Because this algorithm only returns
    one optimal alignment, we will have a systematic way of choosing: high road
    or low road. 
    '''
    maxscore = max(scores)
    if road == 'high':
        if maxscore == scores[3]:
            dirmat[x][y] = 3
        elif maxscore == scores[2]:
            dirmat[x][y] = 2
        elif maxscore == scores[1]:
            dirmat[x][y] = 1
        else:
            dirmat[x][y] = 0
    elif road == 'low':
        if maxscore == scores[1]:
            dirmat[x][y] = 1
        elif maxscore == scores[2]:
            dirmat[x][y] = 2
        elif maxscore == scores[3]:
            dirmat[x][y] = 3
        else:
            dirmat[x][y] = 0
    else:
        # Execution should not reach here.
        raise ValueError('Invalid road type.')
    return maxscore, dirmat

def traceback(seq1, seq2, scoremat, dirmat, start_pos):
    '''
    Find the optimal path through the matrix.
    This function traces a path from the position of the max score to an entry with
    a zero. Each move corresponds to a match, mismatch, or gap in one
    or both of the sequences being aligned. Moves are determined by the direction
    matrix passed to the function.
    WHAT EACH MOVE REPRESENTS
        diagonal: match/mismatch
        up:       gap in sequence 1
        left:     gap in sequence 2
    '''
    # Initialize the sequences
    aligned_seq1 = []
    aligned_seq2 = []
    # Get position of the start, the entry with the highest value
    x, y         = start_pos
    # Use the direction matrix that was filled as we built the score matrix
    move         = dirmat[x][y]
    while move != 0: # If you encounter a zero, end the run. You hit an end-gap.
        if move == 3: # Defined in find_move(), 3 = up
            aligned_seq1.append(seq1[x - 1]) # Get the amino acid from seq 1
            aligned_seq2.append('-') # Skip the amino acid from seq 2
            x -= 1 # Decrement the x position to move to the "previous" seq 1 position
        elif move == 2:
            aligned_seq1.append(seq1[x - 1])
            aligned_seq2.append(seq2[y - 1])
            x -= 1 # Decrement both x and y when there's a match or a mismatch
            y -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[y - 1])
            y -= 1

        move = dirmat[x][y] # Find the next move
    
    # As with most dynamic programming algorithms for sequences, optimal sequences
    # are found in reverse. Reverse them to get the sensible sequence. 
    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))

def alignment_string(aligned_seq1, aligned_seq2):
    '''
    Construct a special string showing identities, gaps, and mismatches.
    This string is printed between the two aligned sequences and shows the
    identities (|), gaps (-), and mismatches (:). As the string is constructed,
    it also counts number of identities, gaps, and mismatches and returns the
    counts along with the alignment string.
    AAGGATGCCTCAAATCGATCT-TTTTCTTGG-
    ::||::::::||:|::::::: |:  :||:|   <-- alignment string
    CTGGTACTTGCAGAGAAGGGGGTA--ATTTGG
    '''
    # Build the string as a list of characters to avoid costly string
    # concatenation.
    idents, gaps, mismatches = 0, 0, 0
    alignment_string = []
    for base1, base2 in zip(aligned_seq1, aligned_seq2):
        if base1 == base2:
            alignment_string.append('|')
            idents += 1
        elif '-' in (base1, base2):
            alignment_string.append(' ')
            gaps += 1
        else:
            alignment_string.append(':')
            mismatches += 1

    return ''.join(alignment_string), idents, gaps, mismatches
