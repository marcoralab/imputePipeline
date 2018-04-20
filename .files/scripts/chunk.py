import sys

def chunker(chrom, chunk, end, start=1, ho=False, exact=False,
            start_offset=False, sliding=False, raw=False, zerobased=False):
    if zerobased and start == 1:
        start = 0
    offset1 = 0
    nc_chunk = chunk
    if start_offset:
        nc_end = end - start
        offset0 = start - 1
    else:
        nc_end = end
        offset0 = -1 if zerobased else 0
    if ho:
        offset1 += 1
    if sliding:
        offset1 += sliding
        chunk -= sliding
        nc_end += sliding*2

    nc = -( -nc_end // nc_chunk )

    rg = lambda x: (chrom,
                    offset0 + chunk * x + 1,
                    offset0 + chunk * (x + 1) + offset1)
    fmt = lambda x: "{}:{}-{}".format(*x)
    ranges = [rg(x) for x in range(nc)]
    if exact:
        if ho:
            end += 1
        ranges[0] = (chrom, start, ranges[0][2])
        ranges[-1] = (*ranges[-1][:2], end)
    ranges = [x for x in ranges if x[2] >= start]
    if raw:
        return ranges
    return [fmt(x) for x in ranges]

def test():
    testchunk = lambda x, chunksize: all([i[2] - i[1] ==
                                         chunksize - 1 for i in x])
    testoverlap = lambda x, overlap: x[0][2] == x[1][1] + overlap - 1
    teststart = lambda x, start: x[0][1] == start
    testend = lambda x, end: x[-1][2] >= end
    testend_exact = lambda x, end: x[-1][2] == end

    # Tests Below

    start = 3
    end = 7117
    chunksize = 3000

    # Test half-closed, zero based, offset, inexact
    chunks = chunker(1, chunksize, end, start, ho=False, exact=False,
                     start_offset=True, sliding=False, zerobased=True,
                     raw=True)
    assert testchunk(chunks, chunksize)
    assert testoverlap(chunks, 0)
    assert teststart(chunks, start)
    assert testend(chunks, end)
    assert not testend_exact(chunks, end)
    chunks = chunker(1, chunksize, chunks[1][2], start, ho=False, exact=False,
                     start_offset=True, sliding=False, zerobased=True,
                     raw=True)
    assert testend(chunks, chunks[1][2])
    assert testend_exact(chunks, chunks[1][2])

    # Test half-open, zero based, offset, inexact
    chunks = chunker(1, chunksize, end, start, ho=True, exact=False,
                     start_offset=True, sliding=False, zerobased=True,
                     raw=True)
    assert testchunk(chunks, chunksize + 1)
    assert testoverlap(chunks, 1)
    assert teststart(chunks, start)
    assert testend(chunks, end + 1)
    assert not testend_exact(chunks, end)

    # Test half-open, zero based, un-offset, inexact
    chunks = chunker(1, chunksize, end, start, ho=True, exact=False,
                  start_offset=False, sliding=False, zerobased=True,
                  raw=True)
    assert testchunk(chunks, chunksize + 1)
    assert testoverlap(chunks, 1)
    assert teststart(chunks, 0)
    assert testend(chunks, end + 1)
    assert not testend_exact(chunks, end)

    # Test half-open, zero based, un-offset, exact
    chunks = chunker(1, chunksize, end, start, ho=True, exact=True,
                     start_offset=False, sliding=False, zerobased=True,
                     raw=True)
    assert testoverlap(chunks, 1)
    assert teststart(chunks, start)
    assert testend(chunks, end + 1)
    assert testend_exact(chunks, end + 1)

    # Test half-closed, 1 based, un-offset, exact
    chunks = chunker(1, chunksize, end, start, ho=False, exact=True,
                     start_offset=False, sliding=False, zerobased=False,
                     raw=True)
    assert testoverlap(chunks, 0)
    assert teststart(chunks, start)
    assert testend(chunks, end)
    assert testend_exact(chunks, end)

    # Test half-closed, 1 based, offset, exact
    chunks = chunker(1, chunksize, end, start, ho=False, exact=True,
                     start_offset=True, sliding=False, zerobased=False,
                     raw=True)
    assert testchunk(chunks[:2], chunksize)
    assert testoverlap(chunks, 0)
    assert teststart(chunks, start)
    assert testend(chunks, end)
    assert testend_exact(chunks, end)

    # Test half-closed, zero based, un-offset, inexact, sliding
    chunks = chunker(1, chunksize, end, start, ho=False, exact=False,
                     start_offset=False, sliding=1000, zerobased=True,
                     raw=True)
    assert testchunk(chunks, chunksize)
    assert testoverlap(chunks, 1000)
    assert teststart(chunks, 0)
    assert testend(chunks, end + 1)
    assert not testend_exact(chunks, end)

    # Test half-open, one based, offset, inexact, sliding
    chunks = chunker(1, chunksize, end, start, ho=True, exact=False,
                     start_offset=True, sliding=1000, zerobased=True,
                     raw=True)
    assert testchunk(chunks, chunksize+1)
    assert testoverlap(chunks, 1001)
    assert teststart(chunks, start)
    assert testend(chunks, end + 2)
    assert not testend_exact(chunks, end)

if __name__ == "__main__":
    chunksize = 20000000
    chrom = str(sys.argv[1])
    start = int(sys.argv[2])
    end = int(sys.argv[3])
    if len(sys.argv) == 5:
        ftype = sys.argv[4]
    else:
        ftype = "human"

    if ftype == "human":
        chunks = chunker(chrom, chunksize, end, start)
        print("\n".join(chunks))
    elif ftype == "bed":
        chunks = chunker(chrom, chunksize, end - 1, start - 1,
                         raw=True, zerobased=True, ho=True)
        [print('\t'.join([str(y) for y in x])) for x in chunks]
    elif ftype == "bcftools":
        import pandas as pd
        chunks = chunker(chrom, chunksize, end, start, raw=True)
        [print('\t'.join([str(y) for y in x])) for x in chunks]
    elif ftype == "bcftools_arg":
        chunks = chunker(chrom, chunksize, end, start)
        print(",".join(chunks))
