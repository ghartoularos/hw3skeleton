'''
I/O for getting pairs from pairs of sequences from files.
'''

def get_pairs(filename):
	with open(filename,'r') as file:
			file = file.readlines()
			pairs = [[seqfile for seqfile in pair.strip().split(' ')] 
					 for pair in file]
			seqpairs = []
			for pair in pairs:
				seqs = []
				for i in range(2):
					lines = open(pair[i],'r').readlines()
					seqs.append(''.join([line.strip() for line in lines 
										 if line[0] != '>']))
				seqpairs.append(seqs)
	return seqpairs
