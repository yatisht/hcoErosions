import argparse
import os
import glob
import anydbm
import string
import subprocess
import sys
import tempfile
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import cStringIO
import pickle
import itertools
import collections

CHAIN_ID_DIR = 'chain_ids/'
CHAIN_HELPERS_DIR = 'chains/' 
GAP_DIR = 'gaps/'
TWO_BIT_TO_FA = 'bin/twoBitToFa'  
TWO_BIT_TO_FA_CDS = 'bin/twoBitToFaCDS'
BIG_BED_TO_BED = 'bin/bigBedToBed'
TWO_BIT_SEQ_DIR = '2bit/'
BED12_TO_BED6 = 'bin/bedToExons'

chain_id_dict = {} 

class GapTrack:
    def __init__(self, bedRegions):
        self.data = collections.defaultdict(list)
        for reg in bedRegions:
            self.data[reg.chrom].append(reg)
        for k in self.data:
            self.data[k].sort()

    def __len__(self):
        return sum([ len(self.data[x]) for x in self.data ])

    def hasGapInRegion(self, bedRegion):
        if bedRegion.chrom not in self.data:
            return False
        curList = self.data[bedRegion.chrom]
        lo = 0
        hi = len(curList) - 1
        while hi - lo > 8:
            mid = (lo + hi) / 2
            if curList[mid].start >= bedRegion.end: hi = mid
            elif curList[mid].end <= bedRegion.start: lo = mid
            else: return True  # This gap must overlap
        for i in range(lo, hi+1):
            if not (curList[i].start >= bedRegion.end or 
                    curList[i].end <= bedRegion.start): return True
        return False

    def numGapsOverlappingRegion(self, bedRegion):
        return len(self.getGapsOverlappingRegion(bedRegion))


    def getGapsOverlappingRegion(self, bedRegion):
        if bedRegion.chrom not in self.data:
            return []
        # Find the rightmost region starting to the left of bedRegion
        curList = self.data[bedRegion.chrom]
        lo = 0  # inclusive index
        hi = len(curList)  # exclusive index
        while hi - lo > 8:
            mid = (lo + hi) / 2
            if curList[mid].start >= bedRegion.start: hi = mid
            else: lo = mid

        # Now iterate rightwards from lo 
        ans = []
        for i in xrange(lo, len(curList)):
            if curList[i].start >= bedRegion.end: 
                break
            if curList[i].end > bedRegion.start: 
                ans.append(curList[i])
        return ans

    def __iter__(self):
        return GapTrackIterator(self)

class GapTrackIterator:
    def __init__(self, gt):
        self.gt = gt
        self.chroms = gt.data.keys()
        self.curChrom = 0
        self.curIndex = 0

    def __iter__(self):
        return self

    def next(self):
        self.curIndex += 1
        if self.curIndex >= len(self.gt.data[self.chroms[self.curChrom]]):
            self.curIndex = 0
            self.curChrom += 1
            if self.curChrom >= len(self.chroms): raise StopIteration
        return self.gt.data[self.chroms[self.curChrom]][self.curIndex]


def readGapTrack(filename):
    return GapTrack(readBedFile(filename))

def readBedFile(filename=None):
    if filename is None:
        f = sys.stdin
    else:
        f = open(filename)
    ans = [ BedRegion.parse(x.strip()) for x in f if x[0] != '#' ]  
            # strip comments
    if filename is not None:
        f.close()
    return ans


class ChainLink:
    def __init__(self, tBed, qBed, frame, chainId, index):
        self.tBed = tBed
        self.qBed = qBed
        self.frame = int(frame)
        self.chainId = int(chainId)
        self.index = int(index)  # Index of the link within the larger chain

    @staticmethod
    def linkKey(chainId, index):
        "Creates a standardized unique identifier for a chain link"
        return '%d_%d' % (chainId, index)

    @staticmethod
    def linkKeyToNums(key):
        return tuple([ int(x) for x in key.split('_') ])

    def toSqlRecord(self):
        return (self.tBed.start, self.qBed.start, self.tBed.length(), 
                self.frame, self.chainId, self.index)

    def toString(self):
        return ' '.join([ str(x) for x in self.toSqlRecord() ])

    @staticmethod
    def parseSqlRecord(record, chain):
        tStart, qStart, length, frame, chainId, index = [ int(x) if x.isdigit()
                else x for x in record ]
        assert chain.id == chainId
        tBed = BedRegion(chain.tChr, tStart, tStart + length, 
                name = ChainLink.linkKey(chainId, index), 
                strand = chain.tStrand)
        qBed = BedRegion(chain.qChr, qStart, qStart + length, 
                name = ChainLink.linkKey(chainId, index), 
                strand = chain.qStrand)
        return ChainLink(tBed, qBed, frame, chainId, index)


    @staticmethod
    def parseString(line, chain):
        "Parse output of toString() back to ChainLink"
        return ChainLink.parseSqlRecord(line.split(), chain)


    def length(self):
        return self.tEnd - self.tStart


class ChainHeader(object):
    def __init__(self, firstLine):
        (c, self.score, self.tChr, self.tSize, self.tStrand, self.tStart,
                self.tEnd, self.qChr, self.qSize, self.qStrand, self.qStart,
                self.qEnd, self.id) = [ int(x) if x.isdigit() else x for x in
                        firstLine.split() ]

        
    def toSqlRecord(self):
        return (self.score, self.tChr, self.tSize, self.tStrand, self.tStart,
                self.tEnd, self.qChr, self.qSize, self.qStrand, self.qStart,
                self.qEnd, self.id)


    def toHeaderString(self):
        "Convert back to header line in such as what is in a .chain file"
        return 'chain ' + ' '.join([ str(x) for x in self.toSqlRecord() ])


    @staticmethod
    def parseSqlRecord(record):
        firstLine = 'chain ' + ' '.join([ str(x) for x in record ])
        return ChainHeader(firstLine)


    def tBed(self):
        "Returns bed region in target covered by this chain"
        return BedRegion(self.tChr, self.tStart, self.tEnd, self.id, 0, 
                self.tStrand)


    def qBed(self):
        "Returns bed region in query covered by this chain"
        return BedRegion(self.qChr, self.qStart, self.qEnd, self.id, 0, 
                self.qStrand)


class BedRegion:
    def __init__(self, chrom, start, end, name=None, score=0, strand=None):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name
        #Score not necessarily an integer. 
        try:
            self.score = int(score)
        except ValueError:
            self.score = 0
        self.strand = strand

    def length(self):
        return self.end - self.start

    def __str__(self):
        ret = '%s\t%d\t%d' % (self.chrom, self.start, self.end)
        if self.strand:
            # If strand, need to fill in name and score
            name = self.name if self.name else '[NO_NAME]'
            ret += '\t%s\t%d\t%s' % (name, self.score, self.strand)
        elif self.name:
            ret += '\t%s' % (self.name)
        return ret

    def __hash__(self):
        return hash(self.toTuple())

    def __eq__(self, other):
        if not isinstance(other, BedRegion):
            return false
        return self.toTuple() == other.toTuple()

    def __lt__(self, other):
        # Sort by chrom, then start, then end
        if self.chrom != other.chrom:
            return self.chrom < other.chrom
        if self.start != other.start:
            return self.start < other.start
        return self.end < other.end

    def toTuple(self):
        return (self.chrom, self.start, self.end, self.name, self.score,
                self.strand)

    def intersect(self, other):
        if self.chrom != other.chrom: return None
        chrom = self.chrom
        start = max(self.start, other.start)
        end = min(self.end, other.end)
        if start > end: return None

        # Use name if only one unique name between the two of them
        names = set([x for x in (self.name, other.name) if x is not None ])
        if len(names) == 1:
            name = names.pop()
        else:
            name = None

        # Use strand if only one unique strand between the two of them
        strands = set([x for x in (self.strand, other.strand) if x is not None])
        if len(strands) == 1:
            strand = strands.pop()
        else:
            strand = None

        return BedRegion(chrom, start, end, name=name, strand=strand)

    @classmethod
    def parse(cls, line):
        vals = line.split()
        if len(vals) < 3:
            raise Exception('Line is not in bed12 format: only %d fields.' %
                    len(vals))
        elif len(vals) == 3:
            (chrom, start, end) = vals
            return cls(chrom, start, end)
        elif len(vals) == 4:
            (chrom, start, end, name) = vals[:4]
            return cls(chrom, start, end, name)
        else:
            (chrom, start, end, name,score,strand) = vals[:6]
            return cls(chrom, start, end, name, score, strand)

    @classmethod
    def parseExons(cls, line):
        words = line.split()
        if (len(words) < 12):
            raise Exception('Line is not in bed12 format: only %d fields.' %
                    len(words)) 
        toRet = []
        chrom = words[0]
        chromStart = int(words[1])
        name = words[3]
        score = words[4]
        strand = words[5]
        blockSizes = words[10].rstrip(',').split(',')
        blockStarts = words[11].rstrip(',').split(',')
        for blockSize, blockStart in itertools.izip(blockSizes, blockStarts):
            exonstart = chromStart + int(blockStart)
            exonend = exonstart + int(blockSize)
            toRet += [cls(chrom, exonstart, exonend, name, score, strand)]
        return toRet

    @classmethod
    def parseRegionList(cls, line):
        try:
            return cls.parseExons(line)
        except Exception:
            return [ cls.parse(line) ]

def reverseAlignedLinks(links):
    nlinks = []
    strandDict = {'-':'+', '+':'-'}
    for link in reversed(links):
	nlinks += [\
	BedRegion(link.chrom, link.start, link.end, link.name,\
	link.score, strandDict[link.strand])
	]
    return nlinks


def makeTemporary(bedRegions):
    "Creates a temporary file whose contents are the given bed regions"
    tmpBed = tempfile.NamedTemporaryFile(suffix='.bed')
    for b in bedRegions:
        print>>tmpBed, str(b)
    tmpBed.flush()
    return tmpBed

def intersectChainAndRegion(reference, qspecies, chainHelper, chainHeader,
    bedregion, verbose=False):

    if chainHeader == None:
	return ([], [])

    if not chainHelper.computeOverhang(chainHeader, bedregion)[0] > 0:
	return ([], [])

    chainLinks = chainHelper.chainLinksOverlappingRegion(bedregion,\
	chainHeader)

    for link in chainLinks:
	assert link.tBed.strand == '+'

    def g(bed, startTrim, endTrim):
	bedlen = bed.end - bed.start
	assert(startTrim + endTrim <= bedlen)
	if bed.strand == '+':
	    nstart = bed.start + startTrim
	    nend = bed.end - endTrim
	else:
	    nstart = bed.start + endTrim
	    nend = bed.end - startTrim
	return BedRegion(bed.chrom,nstart,nend, \
	bed.name, 0, bed.strand)

    lastStartLEbedStart = -1
    for i in reversed(xrange(len(chainLinks))):
	tBed = chainLinks[i].tBed
	if tBed.start <= bedregion.start:
	    lastStartLEbedStart = i
	    break;
    firstEndGEbedEnd = -1
    for i in xrange(len(chainLinks)):
	tBed = chainLinks[i].tBed
	if tBed.end >= bedregion.end:
	    #Mark the end as the next link after this one
	    firstEndGEbedEnd = i+1
	    break;
    if lastStartLEbedStart == -1 or firstEndGEbedEnd == -1:
	return ([], [])

    if verbose:
	print 'Aligning region ', bedregion, lastStartLEbedStart, \
	len(chainLinks), firstEndGEbedEnd

    refAlignedLinks = []
    queryAlignedLinks = []
    for i in xrange(lastStartLEbedStart, firstEndGEbedEnd):
	link = chainLinks[i]
	if verbose:
	    print 'LINK chain ', chainHeader.id, 'ref', link.tBed
	    print 'LINK chain ', chainHeader.id, 'query', link.qBed

	startTrim = 0
	endTrim = 0
	linklen = link.tBed.end - link.tBed.start
	nTbedStart = bedregion.start
	if link.tBed.end <= bedregion.start:
	    refAlignedLinks += [BedRegion(bedregion.chrom,\
	    bedregion.start, bedregion.start, bedregion.name, \
	    bedregion.score, '+')]
	    queryAlignedLinks += [g(link.qBed, linklen, 0)]
	elif link.tBed.start > bedregion.end:
	    refAlignedLinks += [BedRegion(bedregion.chrom,\
	    bedregion.end, bedregion.end, bedregion.name, \
	    bedregion.score, '+')]
	    queryAlignedLinks += [g(link.qBed, 0, linklen)]
	else:
	    endTrim = max(0,link.tBed.end - bedregion.end)
	    startTrim = max(0,bedregion.start - link.tBed.start)
	    refAlignedLinks += [g(link.tBed,startTrim, endTrim)]
	    queryAlignedLinks += [g(link.qBed, startTrim, endTrim)]

    if bedregion.strand == '-':
	refAlignedLinks = reverseAlignedLinks(refAlignedLinks)
	queryAlignedLinks = reverseAlignedLinks(queryAlignedLinks)

    return (refAlignedLinks, queryAlignedLinks)

def getGaplessAlignment_(tSeqs, tRegions, tAssembly,
    qSeqs, qRegions, qAssembly,
	lowQualityAlignment, verbose, edgeBasesToMask):
    assert len(tRegions) == len(qRegions)
    outRef = cStringIO.StringIO()
    outQuery = cStringIO.StringIO()
    indels = []

    for i in xrange(len(tRegions)):
	assert len(tSeqs[i]) == len(qSeqs[i])

	outRef.write(tSeqs[i])

	queryAlignInGap = lowQualityAlignment(qRegions[i])
	if queryAlignInGap:
	    outQuery.write('?' * len(qSeqs[i]))
	else:
	    outQuery.write(qSeqs[i])

	indels += [str(len(tSeqs[i]))]

	isLast = (i == len(tRegions) - 1)
	if not isLast:
	    next = i+1
	    if tRegions[next].strand == '+':
		tinsert_start = tRegions[i].end
		tinsert_end = tRegions[next].start
		region_start = tRegions[0].start
		region_end = tRegions[-1].end
	    else:
		tinsert_start = tRegions[next].end
		tinsert_end = tRegions[i].start
		region_start = tRegions[-1].start
		region_end = tRegions[0].end
	    if qRegions[next].strand == '+':
		qinsert_start = qRegions[i].end
		qinsert_end = qRegions[next].start
	    else:
		qinsert_start = qRegions[next].end
		qinsert_end = qRegions[i].start
	    qinsert = qinsert_end - qinsert_start
	    qdelete = tinsert = tinsert_end - tinsert_start
	    assert qdelete >= 0
	    assert qinsert >= 0

	    qinsert_region = BedRegion(qRegions[i].chrom, qinsert_start,\
	    qinsert_end, qRegions[i].name, 0, '+')
	    tinsert_region = BedRegion(tRegions[i].chrom, tinsert_start,\
	    tinsert_end, tRegions[i].name, 0, '+')

	    insertionInGap = lowQualityAlignment(qinsert_region)

	    if tinsert_end < region_start + edgeBasesToMask:
		insertionInGap = True
	    if tinsert_start >= region_end - edgeBasesToMask:
		insertionInGap = True

	    if insertionInGap:
		outRef.write("-" * qdelete)
		outQuery.write("?" * qdelete)
	    else:
		outRef.write("-" * qdelete)
		outQuery.write("-" * qdelete)

	    if verbose:
		print "gap=", insertionInGap, " when aligning ", tinsert_region,'\n', qinsert_region

	    indels += [getInsDelStr(qinsert, qdelete, insertionInGap)]


    _outRef = list(outRef.getvalue())
    _outQuery = list(outQuery.getvalue())
    outRef.close()
    outQuery.close()

    N = len(_outRef)
    for i in xrange(0, min(N, edgeBasesToMask)):
	_outRef[i] = getMaskedBase(_outRef[i])
	_outQuery[i] = getMaskedBase(_outQuery[i])
    for i in xrange(max(0, N - edgeBasesToMask), N):
	_outRef[i] = getMaskedBase(_outRef[i])
	_outQuery[i] = getMaskedBase(_outQuery[i])

    return ("".join(_outRef), "".join(_outQuery), ",".join(indels))

def getMaskedBase(base):
    if base == '-':
	return base
    return '?'

def getSequence(assembly, bedRegions, twoBitFile=None):
    ret_mapping = {}
    nBedRegions = []
    for i, b in enumerate(bedRegions):
	if b.start == b.end:
	    ret_mapping[i] = -1
	else:
	    ret_mapping[i] = len(nBedRegions)
	    nBedRegions += [b]
    
    if len(nBedRegions) == 0:
	toRet = []
	for i in xrange(len(bedRegions)):
	    toRet += ['']
	return toRet

    tmpBed = makeTemporary(nBedRegions)
    tmpOut = tempfile.NamedTemporaryFile()

    if not twoBitFile:
        twoBitFile = getTwoBitFilename(assembly)
    with open(os.devnull, 'w') as devnull:
        subprocess.check_call([TWO_BIT_TO_FA, twoBitFile, tmpOut.name,
                "-bed=%s" % tmpBed.name, '-noMask'], shell=True)

    seqs = []
    for line in tmpOut:
        line = line.strip()
        if line[0] == '>':
            seqs.append('')
        else:
            seqs[-1] = seqs[-1] + line
    assert len(seqs) == len(nBedRegions)

    toRet = []
    for i in xrange(len(bedRegions)):
	mapping = ret_mapping[i]
	if mapping == -1:
	    toRet += ['']
	else:
	    toRet += [seqs[mapping]]

    # Clean-up
    tmpBed.close()
    tmpOut.close()
    return toRet

def getSequenceBatch(assembly, bedRegions):
    bedRegions_all = []
    for list in bedRegions:
	bedRegions_all += list
    seqs = getSequence(assembly, bedRegions_all)
    toRet = []
    i = 0
    for list in bedRegions:
	toRet += [seqs[i:(i+len(list))]]
	i += len(list)
    return toRet

def getInsDelStr(insertions, deletions, isPoorQuality):
    if isPoorQuality:
	return "?%s-%s" % (str(insertions), str(deletions))
    else:
	return "%s-%s" % (str(insertions), str(deletions))


def getGaplessAlignment(actualRegions, tRegions, tAssembly, qRegions, qAssembly,
	lowQualityAlignment, verbose=False, edgeBasesToMask = 6):
    assert len(tRegions) == len(qRegions)
    
    tSeqs = getSequenceBatch(tAssembly, tRegions)
    qSeqs = getSequenceBatch(qAssembly, qRegions)

    refAlign = []
    queryAlign = []
    indels = []

    for region in xrange(len(actualRegions)):
	if len(tSeqs[region]) == 0:
	    region_len = actualRegions[region].end - actualRegions[region].start
	    #Sequence aligns into gap.
	    refAlign += ["-" * region_len]
	    queryAlign += ["?" * region_len]
	    indels += [getInsDelStr(0,region_len,True)]
	    continue

	(sub_refalign, sub_queryalign, sub_indels) = \
	getGaplessAlignment_(
tSeqs[region], tRegions[region], tAssembly, qSeqs[region],
qRegions[region], qAssembly,
	lowQualityAlignment, verbose, edgeBasesToMask)

	refAlign += [sub_refalign]
	queryAlign += [sub_queryalign]
	indels += [sub_indels]

    return ("".join(refAlign), "".join(queryAlign),\
    '|'.join(indels))

def getRegionAlignmentFromHomologChain(reference, qspecies,
    chainHelper, homologChain, lowQualityAlignment, bedregions, verbose=False):

    refAlignedLinks = []
    queryAlignedLinks = []

    #Use chain helpers to look up homologs of region 
    print bedregions
    for bedregion in bedregions:
	(sub_refAlignedLinks, sub_queryAlignedLinks) = \
	intersectChainAndRegion(reference, qspecies,\
	chainHelper, homologChain, bedregion, verbose)

	refAlignedLinks += [sub_refAlignedLinks]
	queryAlignedLinks += [sub_queryAlignedLinks]


    (refAlignment, queryAlignment, indels) = \
    getGaplessAlignment(bedregions, refAlignedLinks, \
    reference, queryAlignedLinks, qspecies, lowQualityAlignment, verbose)

    return (refAlignment, queryAlignment, indels)


class ChainHelper:
    def __init__(self, chainBigBed, chainDb, linkDb):
        self.chainBigBed = chainBigBed
        self.chainDb = anydbm.open(chainDb)
        self.linkDb = anydbm.open(linkDb)

    @staticmethod
    def initFromDir(dir, refassembly, queryassembly):
        chainBigBed = _findUniqueFile(dir, refassembly + '.' + queryassembly +
                                      '*.chain.bb')
        chainDb = _findUniqueFile(dir, refassembly + '.' + queryassembly +
                                      '*.chain.db')
        linkDb = _findUniqueFile(dir, refassembly + '.' + queryassembly +
                                      '*.link.db')
        return ChainHelper(chainBigBed, chainDb, linkDb)


    def close(self):
        self.chainDb.close()
        self.linkDb.close()


    def getChainHeader(self, chainId):
        print chainId
        return ChainHeader(self.chainDb[chainId])


    def getChainLink(self, chainHeader, index):
        key = ChainLink.linkKey(chainHeader.id, index)
        if key not in self.linkDb: return None
        return ChainLink.parseString(self.linkDb[key], chainHeader)


    def chainsOverlappingRegion(self, bedRegion):
        tmpBed = tempfile.NamedTemporaryFile(suffix='.bed')
        subprocess.check_call([ BIG_BED_TO_BED, self.chainBigBed, tmpBed.name, 
                '-chrom=%s' % bedRegion.chrom, '-start=%d' % bedRegion.start,
                '-end=%d' % bedRegion.end ], shell=True) 
	possibleChains = []
	for line in tmpBed:
	    possibleChains += [BedRegion.parse(line.strip())]
	if len(possibleChains) > 10:
	    print >>sys.stderr, "WARNING: More than 10 chains aligning"\
	    " region ", bedRegion, ". Will ignore all but 10 longest."
	possibleChains = sorted(possibleChains, key=lambda x: -(x.end - x.start))
	chainIds = []
	for i in xrange(min(10,len(possibleChains))):
	    chainIds += [possibleChains[i].name]
        tmpBed.close()
        chainHeaders = [ self.getChainHeader(x) for x in chainIds ]
        return chainHeaders

    def chainLinkBefore(self, bedRegion, chainHeader):
        lo = 1
        hi = 1
        while True:
            hiLink = self.getChainLink(chainHeader, hi)
            if hiLink is None or hiLink.tBed.start > bedRegion.start: break
            hi *= 2
        while hi - lo > 4:
            mid = (lo + hi)/2
            midLink = self.getChainLink(chainHeader, mid)
            if midLink is None or midLink.tBed.start > bedRegion.start: hi = mid
            elif midLink.tBed.start == bedRegion.start: return mid
            else: lo = mid

        # Linear search over the small range remaining
        while True:
            link = self.getChainLink(chainHeader, lo)
            if link is None or link.tBed.start > bedRegion.start: 
                return lo - 1
            lo += 1


    def chainLinksOverlappingRegion(self, bedRegion, chainHeader):
        assert chainHeader.tStrand == '+'
        firstLinkIndex = self.chainLinkBefore(bedRegion, chainHeader)
        linkIndices = [firstLinkIndex]
        lastLinkIndex = firstLinkIndex + 1
        while True:
            lastLink = self.getChainLink(chainHeader, lastLinkIndex)
            if lastLink == None: break
            linkIndices.append(lastLinkIndex)
            if lastLink.tBed.start >= bedRegion.end: break
            lastLinkIndex += 1
        chainLinks = [ self.getChainLink(chainHeader, x) for x in linkIndices ]
        return chainLinks


    @staticmethod
    def computeOverhang(chain, bedRegion):
        leftOverhang = bedRegion.start - chain.tStart
        rightOverhang = chain.tEnd - bedRegion.end
        arr = [ leftOverhang, rightOverhang ]
        return (min(arr), max(arr))  # Sort first by min, then by max


def _findUniqueFile(dir, pattern):
    files = glob.glob(os.path.join(dir, pattern))
    if len(files) == 1: return files[0]
    else:
        raise IOError("Found %d files matching pattern '%s' in directory %s" %
                (len(files), dir, pattern))

def readArgs():
    parser = argparse.ArgumentParser('Computes a diff between a set of'
    ' genes in a reference assembly in a set of other query assemblies')
    parser.add_argument('-reference', '--reference',
            help='Reference assembly name (e.g. hg19)',
	    required=True)
    parser.add_argument('-genes', '--genes',
            help='Bed12 file containing genes in reference',
	    required=True)
    parser.add_argument('-queries','--queries',
            help='Comma separated list of query assembly names',
	    required=True)
    parser.add_argument('-outdir','--outdir',
            help='Directory to output generated files',
	    required=True)
    parser.add_argument('-verbose','--verbose',
            help='Be verbose. Useful for debugging.',
	    action='store_const',
	    const=True,
	    required=False,
	    default=False)
    args = vars(parser.parse_args())
    return args

def getTwoBitFilename(assembly):
    return '%s/%s.2bit' % (TWO_BIT_SEQ_DIR, assembly)

def getCDS(assembly, bedRegions):
    tmpBed = makeTemporary(bedRegions)
    tmpBed6 = tempfile.NamedTemporaryFile()
    tmpOut = tempfile.NamedTemporaryFile()

    twoBitFile = getTwoBitFilename(assembly)

    print bedRegions
    strand = bedRegions[0].split()[5]
    with open(os.devnull, 'w') as devnull:
        subprocess.check_call([BED12_TO_BED6, "-cdsOnly", 
                               "-i", tmpBed.name, "-o", tmpBed6.name], shell=True) 

    with open(os.devnull, 'w') as devnull:
        subprocess.check_call([TWO_BIT_TO_FA, twoBitFile, tmpOut.name,
                "-bed=%s" % tmpBed6.name,'-noMask'], shell=True)
    
    curr_transcript = ''
    curr_exons = ['']
    seqs = []
    for line in tmpOut:
        line = line.strip()
        if line[0] == '>':
            if curr_transcript != line:
                seqs.append('')
                curr_exons = ['']
            elif (curr_transcript == line):
                curr_exons.append('')
            curr_transcript = line
        else:
            curr_exons[-1] = curr_exons[-1] + line
            seq =  ''.join(curr_exons) if (strand == '+') else \
                    ''.join(curr_exons[::-1])
            seqs[-1] = seq
    
    print seqs, bedRegions
    assert len(seqs) == len(bedRegions)

    # Clean-up
    tmpBed.close()
    tmpOut.close()
    tmpBed6.close()
    return seqs

def getChainHelper(refassembly, queryassembly):
    chainHelperDir = '%s' % (CHAIN_HELPERS_DIR)
    toRet = None
    try:
	toRet = ChainHelper.initFromDir(chainHelperDir, refassembly,
                                 queryassembly)
    except:
	print >>sys.stderr, "WARNING: Could not load chains from ",\
	refassembly,' to ', queryassembly
        print  >>sys.stderr, chainHelperDir
    return toRet

def getGapTrackFilename(assembly):
    return '%s/%s.gap.bed' % (GAP_DIR, assembly)

def getGapTrack(assembly, verbose=True):
    if verbose:
	print "Loading gap track for ", assembly
    toRet = readGapTrack(getGapTrackFilename(assembly))
    if verbose:
	print "(Done.)"
    return toRet

def printReferenceGene(reference, line, out, complete_transcripts, verbose=False):
    words = line.split()
    if verbose:
	print line
    cdss = getCDS(reference, [line])
    assert len(cdss) == 1
    if (len(cdss[0]) == 0 or len(cdss[0]) % 3 !=0):
	print >>sys.stderr, 'WARNING: incomplete cds in ', line,\
	    'ignoring transcript.\n'
	return
    cds_aa = str(Seq(cdss[0], generic_dna).translate())
    if (cds_aa[-1] != '*'):
	print >>sys.stderr, 'WARNING: incomplete cds in ', line,\
	    'ignoring transcript.\n'
	return
    out.write('\t'.join([words[3], cds_aa, str(Seq(cdss[0], generic_dna)),\
                         cdss[0],  words[10]]))
    out.write('\n')
    out.flush()

    complete_transcripts += [line]

def _translate(cds):
    return Seq(cds.replace('-','A').replace('?','A')).translate()

def getAADiff(refcds, qcds):
    outAADiff = cStringIO.StringIO()

    assert(len(refcds) == len(qcds))
    assert(len(refcds)%3 == 0)

    refAA = _translate(refcds)
    qAA = _translate(qcds)

    for codon in xrange(len(refAA)):
	i = codon * 3
	#Count mutations in [i, i+3)
	mutations = 0
	isAlignedCodon = True #Does this codon align?
	codonInGap = False #Should we print a '?' for this codon?
	for j in xrange(3):
	    if qcds[i+j] == '-':
		isAlignedCodon = False
	    if qcds[i+j] == '?':
		isAlignedCodon = False
		codonInGap = True
	    if refcds[i+j] != qcds[i+j]:
		mutations += 1
	if isAlignedCodon:
	    if refAA[codon] == qAA[codon]:
		#Synonymous => print AA upper case
		outAADiff.write(qAA[codon].upper())
	    else:
		#Nonsynonymous => print AA lower case
		outAADiff.write(qAA[codon].lower())
	    outAADiff.write(str(mutations))
	else:
	    if codonInGap:
		outAADiff.write('?')
	    else:
		outAADiff.write('-')

    aadiff = outAADiff.getvalue()
    outAADiff.close()

    return aadiff

def chainAlignsRegion(chainHelper, chain, tBed, verbose=False):
    links = chainHelper.chainLinksOverlappingRegion(tBed, chain)
    for link in links:
	if link.tBed.start < tBed.end and \
	    link.tBed.end >= tBed.start:
	    if verbose:
		print 'Chain (competitor?) maps ',tBed
		print 'to ',link.qBed
	    return True
    return False

def getChainsOverlappingGene(reference, qspecies, chainHelper, bedline, verbose=False):
    if chainHelper is None:
	#No chains available.
	return []

    exons = BedRegion.parseExons(bedline)
    if len(exons) == 0:
	print >>sys.stderr, "ERROR: Bedline with no exons, ", bedline
	exit(1)
    generegion = BedRegion.parse(bedline)
    geneLen = generegion.end - generegion.start
    transcript_id = bedline.split()[3]
    print 'Transcript:', transcript_id

    #This only returns the top 10 longest chains overlapping this
    #region, but that's fine
    overlappingChainIds = chain_id_dict.get(transcript_id, [])

    S = []
    P = []
    for chainId in overlappingChainIds:
	chain = chainHelper.getChainHeader(chainId)
	#Calculate chain overlap of gene in reference
	#and filter out chains with short overlap
	chainOverlapStart = min(
	max(generegion.start, chain.tStart),generegion.end)
	chainOverlapEnd = min(
	max(generegion.start, chain.tEnd),generegion.end)
	chainOverlapLen = chainOverlapEnd - chainOverlapStart

	#Calculate span of the chain in the query
	#This is an upper bound on the length of the aligned-to region
	#according to this chain.

	chainQuerySpan = chain.qEnd - chain.qStart

	S += [(chainId, chain)]

	#Compute how many exons this chain aligns at least partially
	numOverlappingBases = 0
	numOverlappingExons = 0
	for exon in exons:
	    exon_len = exon.end - exon.start
	    if chainAlignsRegion(chainHelper, chain, exon):
		numOverlappingBases += exon_len
		numOverlappingExons += 1

	P += [(chainId, chain)]

    toRet = None

    #Do we have any chains?
    if len(S) == 0:
	pass
    #Otherwise, is P empty?
    elif len(P) == 0:
	toRet = S[0]
    #Or, do we have multiple possible orthologs? (confusion due to paralogs)
    elif len(P) > 1:
	pass
    #Remaining case is we have one good ortholog
    else:
	toRet = P[0]

    if toRet is None:
	return []

    if verbose:
	print 'Ortholog chain being used:', toRet[0], toRet[1].tChrom, \
	toRet[1].tStart, toRet[1].tEnd

    #Return just the chain (not the id)
    return [toRet[1]]

def getPairwiseGeneDiffs(reference, qspecies, chainHelper, queryGapTrack,
    bedline, verbose=False):
    bedWords = bedline.split()
    cdsStart = int(bedWords[6])
    cdsEnd = int(bedWords[7])
    strand = bedWords[5]
    name = bedWords[3]
    geneRegion = BedRegion.parse(bedline)

    overlappingChains = getChainsOverlappingGene(\
	reference, qspecies, chainHelper, \
	bedline, verbose)

    if len(overlappingChains) == 0:
	orthologChain = None
    else:
	orthologChain = overlappingChains[-1]

    def lowQualityAlignment(qBed):
	if qBed is not None:
	    if queryGapTrack.hasGapInRegion(qBed):
		return True
	return False


    exons = BedRegion.parseExons(bedline)
    if strand == '-':
	exons = [a for a in reversed(exons)]

    exon_cds = []

    for exon in exons:
	nstart = min(max(exon.start,cdsStart),exon.end)
	nend = max(min(exon.end,cdsEnd), exon.start)
	exon_cds_len = nend - nstart
	assert exon_cds_len >= 0
	if exon_cds_len == 0:
	    continue
	exon_cds += [BedRegion(exon.chrom, nstart, nend, name,
	exon.score, exon.strand)]

    return getRegionAlignmentFromHomologChain(reference, \
    qspecies, chainHelper, orthologChain, lowQualityAlignment, exon_cds, verbose)

def printPairwiseGeneDiffs(reference, qspecies, chainHelper,
queryGapTrack, bedline, out, verbose=False):
    gene = BedRegion.parse(bedline)


    (refAlignment, queryAlignment, indels) = \
    getPairwiseGeneDiffs(reference, qspecies, chainHelper, queryGapTrack, bedline, verbose)

    if len(refAlignment) == 0:
	return

    if len(refAlignment)%3 != 0:
	print >>sys.stderr, 'Refcds not length multiple of 3',\
	    bedline
	exit(1)

    aadiff = getAADiff(refAlignment, queryAlignment)

    out.write("%s\t%s\t%s\t" % (qspecies, gene.name, indels))
    out.write("%s\t%s\n" % (aadiff, queryAlignment))
    out.flush()

def printGeneDiffs(reference, genes, queries, outdir, verbose):
    complete_transcripts = []

    with open(os.path.join(outdir,\
    ''.join([reference, '.', 'referenceGenes.out'])), 'w')\
    as refgenefile:
	with open(genes, "r") as genes_:
	    for line in genes_:
		printReferenceGene(reference, line, refgenefile,
		complete_transcripts, verbose);
    for qassembly in queries.rstrip(',').split(","):
	queryGapTrack = getGapTrack(qassembly)
	chainHelper = getChainHelper(reference, qassembly)
	with open(os.path.join(outdir,\
	''.join([reference,'.',qassembly,'.codingMutations.out'])), 'w')\
	as difffile:
	    for line in complete_transcripts:
		print "Processing", line, "query=",qassembly
		printPairwiseGeneDiffs(reference, qassembly,\
		chainHelper, queryGapTrack, line, difffile, verbose);

def tracefunc(frame, event, arg, indent=[0]):
    if event == "call":
        indent[0] += 2
        print("-" * indent[0] + "> call function", frame.f_code.co_name)
    elif event == "return":
        print("<" + "-" * indent[0], "exit function",
              frame.f_code.co_name)
        indent[0] -= 2
    return tracefunc

if __name__ == "__main__":
    args = readArgs() 
    chain_id_dict = pickle.load(open(CHAIN_ID_DIR + args['queries'] + '.p')) 
    printGeneDiffs(**args)
