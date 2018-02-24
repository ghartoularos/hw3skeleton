from hw3skeleton import smithwaterman, pairsio


def test_smith_waterman():
	seq1 = 'PTKLAVIGAGAVGSTLAFAAAQRGIAREIVLEDIAKE' + \
	'RVEAEVLDMQHGSSFYPTVSIDGSDDPEICRDADMVV'

	seq2 = 'ITAGPRQKPGQSRLELVGATVNILKAIMPNLVKVAPN' + \
	'AIYMLITNPVDIATHVAQKLTGLPENQIFGSG'

	smithwaterman.main('submats/BLOSUM50', seq1, seq2)
	return