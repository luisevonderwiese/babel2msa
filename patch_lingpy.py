from lingpy.algorithm import * 

def my_corrdist(
        threshold,
        seqs,
        gops,
        pros,
        gop,
        scale,
        factor,
        scorer,
        mode,
        restricted_chars
        ):
    """
    Create a correspondence distribution for a given language pair.

    Parameters
    ----------
    threshold : float
        The threshold of sequence distance which determines whether a sequence
        pair is included or excluded from the calculation of the distribution.
    seqs : list
        The sequences passed as a two-dimensional of sequence pairs.
    gops : list
        The gap opening penalties, passed as individual lists of penalties for each
        sequence.
    pros : list
        The of prosodic strings for each sequence.
    gop : int
        The general gap opening penalty which will be used to introduce a gap
        between the two profiles.
    scale : float
        The gap extension scale by which consecutive gaps are reduced. LingPy
        uses a scale rather than a constant gap extension penalty.
    factor : float
        The factor by which matches are increased when two segments occur in
        the same prosodic position of an alignment.
    scorer : { dict, :py:class:`lingpy.algorithm.cython.misc.ScoreDict` }
        The scoring function which needs to provide scores for all
        segments in the two profiles.
    mode : { "global", "local", "overlap", "dialign" }
        Select one of the four basic modes for alignment analyses.
    restricted_chars : { }
        The string containing restricted characters. Restricted characters
        occur, as a rule, in the prosodic strings, not in the normal sequence.
        They need to be computed by computing a consensus string from all
        prosodic strings in the profile.

    Notes
    -----
    This function is the core of the
    :py:class:`~lingpy.compare.lexstat.LexStat` function to compute
    distributions of aligned segment pairs.

    Returns
    -------
    results : tuple
        A dictionary containing the distribution, and the number of included
        sequences.
    """

    # basic defs
# [autouncomment]     cdef int i,j,M,N,lP,l
# [autouncomment]     cdef list seqA,seqB,almA,almB
# [autouncomment]     cdef float sim
    corrs = {}

    # return number of sequences considered for initial distribution
    included = 0

    # get basic params
    lP = len(seqs)

    # check for restricted prostrings

    # carry out alignments
    for i in range(lP):
        # get sequences
        seqA,seqB = seqs[i][0],seqs[i][1]

        # get length of seqs
        M,N = len(seqA),len(seqB)

        # get gops
        gopA = [gop * gops[i][0][j] for j in range(M)]
        gopB = [gop * gops[i][1][j] for j in range(N)]

        # get pros
        proA,proB = pros[i][0],pros[i][1]

        # check for restricted chars
        if not set(restricted_chars).intersection(proA+proB):
            if mode == "global":
                almA,almB,sim = calign.globalign(
                       seqA,
                       seqB,
                       gopA,
                       gopB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer
                       )
            elif mode == "local":
                almA,almB,sim = calign.localign(
                       seqA,
                       seqB,
                       gopA,
                       gopB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer
                       )
                almA = almA[1]
                almB = almB[1]


            elif mode == "overlap":
                almA,almB,sim = calign.semi_globalign(
                       seqA,
                       seqB,
                       gopA,
                       gopB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer,
                       )

            elif mode == "dialign":
                almA,almB,sim = calign.dialign(
                       seqA,
                       seqB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer,
                       )

        else:
            if mode == "global":
                almA,almB,sim = calign.secondary_globalign(
                       seqA,
                       seqB,
                       gopA,
                       gopB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer,
                       restricted_chars
                       )
            elif mode == "local":
                almA,almB,sim = calign.secondary_localign(
                       seqA,
                       seqB,
                       gopA,
                       gopB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer,
                       restricted_chars
                       )
                almA = almA[1]
                almB = almB[1]

            elif mode == "overlap":
                almA,almB,sim = calign.secondary_semi_globalign(
                       seqA,
                       seqB,
                       gopA,
                       gopB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer,
                       restricted_chars
                       )

            elif mode == "dialign":
                almA,almB,sim = calign.secondary_dialign(
                       seqA,
                       seqB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer,
                       restricted_chars
                       )

        # calculate distances
        simA = sum([(1.0 + factor) * scorer[seqA[i],seqA[i]] for i in range(M)])
        simB = sum([(1.0 + factor) * scorer[seqB[i],seqB[i]] for i in range(N)])
        if simA + simB == 0:
            dist = 0.0
        else:
            dist = 1 - ( ( 2 * sim ) / ( simA + simB ) )

        if dist <= threshold:
            included += 1
            l = len(almA)
            for j in range(l):
                try:
                    corrs[almA[j],almB[j]] += 1
                except:
                    corrs[almA[j],almB[j]] = 1

    return corrs, included


calign.corrdist = my_corrdist
