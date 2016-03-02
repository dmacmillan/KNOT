import argparse, sys, os, math
import cPickle as pickle
from copy import deepcopy
#from pybedtools import *
import pysam

def kleat_int(thing):
    try:
        return int(thing)
    except ValueError:
        return None

class BedLine:
    
    def __init__(self, chrom, start, stop, name=None, score=None, strand=None, thickStart=None, thickEnd=None, itemRgb=None, blockCount=None, blockSizes=None, blockStarts=None):
        self.chrom = chrom
        self.chromStart = int(start)
        self.chromEnd = int(stop)
        self.name = name
        self.score = score
        self.strand = strand
        self.thickStart = thickStart
        self.thickEnd = thickEnd
        self.itemRgb = itemRgb
        self.blockCount = blockCount
        self.blockSizes = blockSizes
        self.blockStarts = blockStarts
        if score:
            self.score = int(score)
        self.strand = strand
        if thickStart:
            self.thickStart = int(thickStart)
        if thickEnd:
            self.thickEnd = int(thickEnd)
        self.itemRgb = itemRgb
        if blockCount:
            self.blockCount = int(blockCount)
        if blockSizes:
            self.blockSizes = [int(x) for x in filter(None, blockSizes.split(','))]
        if blockStarts:
            self.blockStarts = [int(x) for x in filter(None, blockStarts.split(','))]

    def __str__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}'.format(self.chrom, self.chromStart, self.chromEnd, self.name, self.score, self.strand)

class KleatResult:

    def __init__(self, gene, transcript, transcript_strand, coding, contig, chromosome, cleavage_site, within_UTR, distance_from_annotated_site, ESTs, length_of_tail_in_contig, number_of_tail_reads, number_of_bridge_reads, max_bridge_read_tail_length, bridge_read_identities, tail_and_bridge_reads, number_of_link_pairs, max_link_pair_length, link_pair_identities, pas, utr3):
        self.gene = gene
        self.transcript = transcript
        self.transcript_strand = transcript_strand
        self.coding = coding
        self.contig = contig
        self.chromosome = chromosome
        self.cleavage_site = int(cleavage_site)
        if (within_UTR == 'no'):
            self.within_UTR = False
        else:
            self.within_UTR = True
        self.distance_from_annotated_site = kleat_int(distance_from_annotated_site)
        self.ESTs = ESTs
        self.length_of_tail_in_contig = kleat_int(length_of_tail_in_contig)
        self.number_of_tail_reads = kleat_int(number_of_tail_reads)
        self.number_of_bridge_reads = kleat_int(number_of_bridge_reads)
        self.max_bridge_read_tail_length = kleat_int(max_bridge_read_tail_length)
        self.bridge_read_identities = bridge_read_identities
        self.tail_and_bridge_reads = kleat_int(tail_and_bridge_reads)
        self.number_of_link_pairs = kleat_int(number_of_link_pairs)
        self.max_link_pair_length = kleat_int(max_link_pair_length)
        self.link_pair_identities = link_pair_identities
        self.pas = self.utr3 = None
        if (pas != '-'):
            self.pas = [int(x) for x in pas.split(':')]
        if (utr3 != '-'):
            self.utr3 = [int(x) for x in utr3.split('-')]

    def __str__(self):
        atts = [self.gene, self.transcript, self.transcript_strand, self.coding, self.contig, self.chromosome, self.cleavage_site, self.within_UTR, self.distance_from_annotated_site, self.ESTs, self.length_of_tail_in_contig, self.number_of_tail_reads, self.number_of_bridge_reads, self.max_bridge_read_tail_length, self.bridge_read_identities, self.tail_and_bridge_reads, self.number_of_link_pairs, self.max_link_pair_length, self.link_pair_identities, self.pas, self.utr3]
        atts = [str(x) for x in atts]
        return ('\t').join(atts)

def parseKleat(kleat, min_bridge_read_tail_len=None, min_num_bridge_reads=None, min_tail_len=None, min_num_tail_reads=None, with_pas=False):
    results = []
    with open(kleat, 'r') as f:
        f.readline()
        for line in f:
            result = KleatResult(*line.strip().split('\t'))
            if min_bridge_read_tail_len and (result.max_bridge_read_tail_length < min_bridge_read_tail_len) and not result.number_of_tail_reads:
                continue
            elif min_num_bridge_reads and (result.number_of_bridge_reads < min_num_bridge_reads) and not result.number_of_tail_reads:
                continue
            elif with_pas and not result.pas:
                continue
            results.append(result)
    return results

class Region:
    def __init__(self, start=None, end=None, piles=None, length=None):
        try:
            self.start = int(start)
        except (ValueError, TypeError) as e:
            self.start = None
        try:
            self.end = int(end)
        except (ValueError, TypeError) as e:
            self.end = None
        self.piles = piles
        self.length = length

class GTF:
    def __init__(self, seqname=None, source=None, feature=None, start=None, end=None, score=None, strand=None, frame=None, attribute=None):
        self.seqname = seqname
        self.source = source
        self.feature = feature
        try:
            self.start = int(start)
        except (ValueError, TypeError) as e:
            self.start = None
        try:
            self.end = int(end)
        except (ValueError, TypeError) as e:
            self.end = None
        try:
            self.score = int(score)
        except (ValueError, TypeError) as e:
            self.score = None
        self.strand = strand
        self.frame = frame
        self.attribute = attribute

    def __str__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(self.seqname, self.source, self.feature, self.start, self.end, self.score, self.strand, self.frame, self.attribute)

def parseGTF(gtffile, seqnames=None, sources=None, features=None , not_sources=None, not_features=None, gzipped=False, add_chr=True):
    results = []
    if gzipped:
        f = gzip.open(gtffile, 'rb')
    else:
        f = open(gtffile, 'r')
    for line in f:
        if line[0] == '#':
            continue
        gtf = GTF(*line.strip().split('\t'))
        if add_chr:
            gtf.seqname = 'chr' + gtf.seqname
        if seqnames and gtf.seqname not in seqnames:
            continue
        attributes = {}
        for attr in [x.split() for x in gtf.attribute.split(';')][:-1]:
            attributes[attr[0]] = attr[1][1:-1]
        gtf.attribute = attributes
        if not_sources and gtf.source in not_sources:
            continue
        elif not_features and gtf.feature in not_features:
            continue
        if sources and gtf.source not in sources:
            continue
        elif features and gtf.feature not in features:
            continue
        results.append(gtf)
    f.close()
    return results

def mergeGTFList(gtf_list):
    res = [gtf_list[0]]
    for i in xrange(1,len(gtf_list)):
        if (gtf_list[i].start <= res[-1].end):
            res[-1] = mergeTwoGTFs(gtf_list[i],res[-1])
        else:
            res.append(gtf_list[i])
    return res

def mergeTwoGTFs(g1, g2, delim='|'):
    res = GTF()
    res.seqname = g1.seqname
    res.source = g1.source + delim + g2.source
    res.feature = g1.feature + delim + g2.feature
    res.start = min(g1.start, g2.start)
    res.end = max(g1.end, g2.end)
    try:
        res.score = float(g1.score) + g2.score / 2
    except TypeError:
        res.score = None
    res.strand = g1.strand
    res.frame = g1.frame
    res.attribute = g1.attribute
    return res

def groupPysamGTF(gtf):
    result = {}
    for g in gtf:
        if g.seqname not in result:
            result[g.seqname] = {g.attribute['gene_name']: [g]}
        if g.attribute['gene_name'] not in result[g.seqname]:
            result[g.seqname][g.attribute['gene_name']] = [g]
        else:
            result[g.seqname][g.attribute['gene_name']].append(g)
    for chrom in result:
        for gene in result[chrom]:
            result[chrom][gene].sort(key=lambda x: x.start)
            result[chrom][gene] = mergeGTFList(result[chrom][gene])
    return result

def groupGTF(gtf):
    result = {}
    for g in gtf:
        if g.seqname not in result:
            result[g.seqname] = {g.attribute['gene_name']: [g]}
        if g.attribute['gene_name'] not in result[g.seqname]:
            result[g.seqname][g.attribute['gene_name']] = [g]
        else:
            result[g.seqname][g.attribute['gene_name']].append(g)
    for chrom in result:
        for gene in result[chrom]:
            result[chrom][gene].sort(key=lambda x: x.start)
            result[chrom][gene] = mergeGTFList(result[chrom][gene])
    return result

def groupKleat(parsed):
    results = {}
    for r in parsed:
        if r.chromosome not in results:
            results[r.chromosome] = {r.gene: [r]}
        if r.gene not in results[r.chromosome]:
            results[r.chromosome][r.gene] = [r]
        else:
            results[r.chromosome][r.gene].append(r)
    return results

def parseConfig(config):
    with open(config, 'r') as f:
        lines = f.readlines()
        lines = [x.strip() for x in lines if x[0] != '#']
    return lines

def mergeKleatResults(sites):
    d = {'cleavage_site': 0,
         'max_len_br': 0,
         'num_br': 0,
         'len_tail_contig': 0,
         'num_tr': 0}
    count = 0
    for c in sites:
        count += 1
        if c.max_bridge_read_tail_length > d['max_len_br']:
            d['max_len_br'] = c.max_bridge_read_tail_length
        d['num_br'] += c.number_of_bridge_reads
        d['cleavage_site'] += c.cleavage_site
        if c.length_of_tail_in_contig > d['len_tail_contig']:
            d['len_tail_contig'] = c.length_of_tail_in_contig
        d['num_tr'] += c.number_of_tail_reads
    d['cleavage_site'] /= count
    res = deepcopy(sites[0])
    res.cleavage_site = d['cleavage_site']
    res.length_of_tail_in_contig = d['len_tail_contig']
    res.number_of_bridge_reads = d['num_br']
    res.number_of_tail_reads = d['num_tr']
    res.max_bridge_read_tail_length = d['max_len_br']
    return res

def kleatLinkage(sites, window=20):
    length = len(sites)
    if length > 1:
        _min = float('inf')
        r = s = None
        for i in xrange(length - 1):
            for j in xrange(i + 1, length):
                dist = abs(sites[i].cleavage_site - sites[j].cleavage_site)
                if dist < _min:
                    r = i
                    s = j
                    _min = dist
        #sites[r] = _min
        if _min <= window:
            sites[r] = mergeKleatResults([sites[r],sites[s]])
            del(sites[s])
            kleatLinkage(sites)
    return sites

def genTrackLine(name, description=None, _type=None, visibility=2, color=None):
    result = ['name="{}"'.format(name)]
    if description:
        result.append('description="{}"'.format(description))
    if _type:
        result.append('type="{}"'.format(_type))
    if visibility:
        result.append('visibility="{}"'.format(visibility))
    if color:
        result.append('color="{}"'.format(color))
    return 'track ' + (' ').join(result) + '\n'

def allDistances(_list):
    result = []
    lenl = len(_list)
    for i in xrange(lenl-1):
        for j in xrange(i+1, lenl):
            dist = abs(_list[i] - _list[j])
            result.append(dist)
    return result

def mean(data):
    """Return the sample arithmetic mean of data."""
    n = len(data)
    if n < 1:
        raise ValueError('mean requires at least one data point')
    return sum(data)/n # in Python 2 use sum(data)/float(n)

def _ss(data):
    """Return sum of square deviations of sequence data."""
    c = mean(data)
    ss = sum((x-c)**2 for x in data)
    return ss

def sstdev(data):
    """Calculates the population standard deviation."""
    n = len(data)
    if n < 2:
        raise ValueError('variance requires at least two data points')
    ss = _ss(data)
    pvar = ss/(n-1) # the population variance
    return pvar**0.5

def pstdev(data):
    """Calculates the population standard deviation."""
    n = len(data)
    if n < 2:
        raise ValueError('variance requires at least two data points')
    ss = _ss(data)
    pvar = ss/n # the population variance
    return pvar**0.5

def writeFile(path, name=None, *lines):
    if name:
        result = os.path.join(path,name)
    else:
        result = path
    with open(result, 'w') as f:
        f.write(('\n').join(lines))
    return result

def computeRatios2(results, annot):
	ratios = {}
	for chrom in results:
		for gene in results[chrom]:
			if gene not in ratios:
				ratios[gene] = {}
			keys = results[chrom][gene].keys()
			strand = annot[chrom][gene][0].strand
			if strand == '+':
				keys = sorted(keys, key=lambda x: int(x.split('-')[0]))
			else:
				keys = sorted(keys, key=lambda x: int(x.split('-')[0]), reverse=True)
			for span in keys:
				for align in results[chrom][gene][span]:
					if align not in ratios[gene]:
						ratios[gene][align] = [results[chrom][gene][span][align]['med']]
					else:
						ratios[gene][align].append(results[chrom][gene][span][align]['med'])
	return ratios

def analyzeRatios(ratios):
    results = {}
    for gene in ratios:
        for sample in ratios[gene]:
            num_regions = len(ratios[gene][sample])
            if num_regions not in results:
                results[num_regions] = 0.5
            else:
                results[num_regions] += 0.5
    return results

def computeRatios(dic, alns, annot):
    for a in alns:
        for chrom in dic:
            alns[a]['changes'][chrom] = {}
            for gene in dic[chrom]:
                #########################################
                # REMOVE THE ENCLOSED CODE FOR ACTUAL RUN
                if len(dic[chrom][gene]) > 2:
                    continue
                #########################################
                alns[a]['changes'][chrom][gene] = []
                strand = annot[chrom][gene][0].strand
                keys = dic[chrom][gene].keys()
                if strand == '+':
                    keys = sorted(keys, key = lambda x: int(x.split('-')[0]))
                else:
                    keys = sorted(keys, key = lambda x: int(x.split('-')[0]), reverse=True)
                for span in keys:
                    try:
                        alns[a]['changes'][chrom][gene].append(dic[chrom][gene][span][a]['med'])
                    except KeyError as e:
                        print chrom, gene, span, a
                        print e
                changes = alns[a]['changes'][chrom][gene]
                dists = []
                for i in xrange(len(changes)-1):
                    try:
                        dist = changes[i+1]/changes[i]
                    except ZeroDivisionError:
                        #dist = float('inf')
                        continue
                    #dist = abs(changes[i] - changes[i+1])
                    dists.append(dist)
                total = sum(dists)

def genResults(annot, kleats):
    results = {}
    fasta = regions = ''
    for chrom in annot:
        if chrom not in kleats:
            continue
        results[chrom] = {}
        for gene in annot[chrom]:
            if gene not in kleats[chrom]:
                continue
            if (not annot[chrom][gene]):
                continue
            results[chrom][gene] = {}
            strand = annot[chrom][gene][0].strand
            gene_start = annot[chrom][gene][0].start
            gene_end = annot[chrom][gene][-1].end
            for a in aligns:
                read_count = 0
                for read in aligns[a]['align'].fetch(chrom, gene_start, gene_end):
                    if (gene_start <= read.pos <= gene_end):
                        read_count += 1
                aligns[a]['read_count'] = read_count
            if strand == '-':
                annot[chrom][gene] = annot[chrom][gene][:1]
            else:
                annot[chrom][gene] = annot[chrom][gene][-2:]
            for region in annot[chrom][gene]:
                last = region.start
                intervals = []
                splices = 0
                cleaved = False
                for k in kleats[chrom][gene]:
                    if (region.start < k.cleavage_site < region.end):
                        cleaved = True
                if not cleaved:
                    continue
                temp = []
                for k in kleats[chrom][gene]:
                    cs = k.cleavage_site
                    key = '{}:{}-{}'.format(chrom,last,cs)
                    if (cs > region.end):
                        cs = region.end
                        pass
                    if cs - last < 20:
                        continue
                    regions += ('\t').join([chrom, str(last), str(cs), '{}_utr_{}'.format(gene,splices), '0', region.strand, '\n'])
                    header = '>' + ('|').join([chrom, gene, str(last), str(cs), region.strand]) + '\n'
                    seq = ref.fetch(chrom, last, cs).upper() + '\n'
                    fasta += header + seq
                    span = '{}-{}'.format(last,cs)
                    results[chrom][gene][span] = {}
                    for a in aligns:
                        m = []
                        for p in aligns[a]['align'].pileup(chrom, last, cs):
                            if (last <= p.pos <= cs):
                                m.append(p.nsegments)
                        sm = sorted(m)
                        lenm = len(sm)
                        if lenm == 0:
                            continue
                        _min = sm[0]
                        _max = sm[-1]
                        med = calcMedian(sm)
                        q1, q3 = calcIQR(sm)
                        _mean = sum(sm)/lenm
                        se = sstdev(sm)/math.sqrt(lenm)
                        temp.append([gene,strand,a,key,lenm,_min,q1,med,q3,_max,_mean,se])
                        #med = sum(m)/lenm
                        if (aligns[a]['read_count'] == 0):
                            print 'No reads: {}'.format(gene)
                            continue
                        #med = float(med)/aligns[a]['read_count']
                        results[chrom][gene][span][a] = {'med': med,
                                                         'read_count': aligns[a]['read_count'],
                                                         'meds': m}
                    if strand == '+':
                        temp = sorted(temp, key=lambda x: int(x[3].split(':')[1].split('-')[0]))
                    else:
                        temp = sorted(temp, key=lambda x: int(x[3].split(':')[1].split('-')[0]), reverse=True)
                    for i in xrange(len(temp)):
                        temp[i] = ('\t').join([str(x) for x in temp[i]])
                        data.append(temp[i])
                    last = cs+1
                    splices += 1
    return results, fasta, regions, data

def calcMedian(lst):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = (lstLen - 1) // 2
    if (lstLen % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1])/2.0

def calcIQR(array):
    length = len(array)
    midpoint = length/2
    if (length % 2 == 0):
        lowermed = calcMedian(array[:midpoint])
        uppermed = calcMedian(array[midpoint:])
    else:
        lowermed = calcMedian(array[:midpoint])
        uppermed = calcMedian(array[midpoint+1:])
    return [lowermed, uppermed]

def computeStats(annot, kleats, aligns):
    data = []
    for chrom in annot:
        if chrom not in kleats:
            continue
        for gene in annot[chrom]:
            if gene not in kleats[chrom]:
                continue
            if (not annot[chrom][gene]):
                continue
            strand = annot[chrom][gene][0].strand
            gene_start = annot[chrom][gene][0].start
            gene_end = annot[chrom][gene][-1].end
            if strand == '-':
                annot[chrom][gene] = annot[chrom][gene][:1]
            else:
                annot[chrom][gene] = annot[chrom][gene][-2:]
            for region in annot[chrom][gene]:
                last = region.start
                intervals = []
                splices = 0
                cleaved = False
                for k in kleats[chrom][gene]:
                    if (region.start < k.cleavage_site < region.end):
                        cleaved = True
                if not cleaved:
                    continue
                temp = []
                for k in kleats[chrom][gene]:
                    cs = k.cleavage_site
                    key = '{}:{}-{}'.format(chrom,last,cs)
                    _len = cs-last
                    if (cs > region.end):
                        cs = region.end
                        pass
                    if cs - last < 20:
                        continue
                    for a in aligns:
                        m = []
                        for p in aligns[a]['align'].pileup(chrom, last, cs):
                            if (last <= p.pos <= cs):
                                m.append(p.nsegments)
                        sm = sorted(m)
                        lenm = len(sm)
                        if lenm == 0:
                            continue
                        _min = sm[0]
                        _max = sm[-1]
                        med = calcMedian(sm)
                        q1, q3 = calcIQR(sm)
                        _mean = sum(sm)/lenm
                        se = sstdev(sm)/math.sqrt(lenm)
                        #data.append(('\t').join([str(x) for x in [gene,strand,a,key,lenm,_min,q1,med,q3,_max,_mean,se]]))
                        temp.append([gene,strand,a,key,lenm,_min,q1,med,q3,_max,_mean,se])
                    last = cs+1
                if strand == '+':
                    temp = sorted(temp, key=lambda x: int(x[3].split(':')[1].split('-')[0]))
                else:
                    temp = sorted(temp, key=lambda x: int(x[3].split(':')[1].split('-')[0]), reverse=True)
                for i in xrange(len(temp)):
                    temp[i] = ('\t').join([str(x) for x in temp[i]])
                    data.append(temp[i])
    return data

# Begin main thread
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Compare two bedgraph tracks")
    
    parser.add_argument('annotation', help='Genome annotation file if GTF format')
    parser.add_argument('kleats', help='File containing list of kleat output files to use')
    parser.add_argument('alignments', help='File containing list of alignment output files to use. Must be in BAM or SAM format. Names of each BAM/SAM file must also differ')
    parser.add_argument('-cw', '--cluster_window', type=int, default=20, help='Set the window size for clustering KLEAT cleavage sites. Default = 20')
    parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Directory to output to. Default is current directory')
    parser.add_argument('-r', '--reference', default='/home/dmacmillan/references/hg19/hg19.fa', help='Path to the reference genome from which to fetch sequences')
    
    args = parser.parse_args()

    ref = pysam.FastaFile(args.reference)

    all_kleats = []

    aligns = {}

    for b in parseConfig(args.alignments):
        name = os.path.basename(os.path.abspath(b)).split('.')[0]
        aligns[name] = {'align': pysam.AlignmentFile(b),
                        'read_count': None,
                        'med': None,
                        'changes': {}}

    for k in parseConfig(args.kleats):
        kleat = parseKleat(k, with_pas=True, min_num_bridge_reads=2, min_bridge_read_tail_len=4)
        all_kleats += kleat

    all_kleats = sorted(all_kleats, key=lambda x: x.cleavage_site)

    kleats = groupKleat(all_kleats)

    for chrom in kleats:
        for gene in kleats[chrom]:
            sites = kleats[chrom][gene]
            sites = kleatLinkage(sites, args.cluster_window)

    #sys.stdout.write('Parsing GTF...')
    #sys.stdout.flush()
    annot = parseGTF(args.annotation, seqnames=['chr{}'.format(x) for x in range(1,23)] + ['chrX', 'chrY'], sources=['protein_coding'], features='UTR')
    #print 'DONE'

    #sys.stdout.write('Grouping GTF...')
    #sys.stdout.flush()
    annot = groupGTF(annot)
    #print 'DONE'

    #data = computeStats(annot, kleats, aligns)
    #temp.append([gene,strand,a,key,lenm,_min,q1,med,q3,_max,_mean,se])
    #print ('\t').join(['GENE','STRAND','SAMPLE','REGION','LENGTH','MIN','Q1','MED','Q3','MAX','MEAN','SE'])
    #print ('\n').join(data)

    regions = ''
    fasta = ''

    results = {}

    sys.stdout.write('Computing coverage of regions...')
    sys.stdout.flush()

    saved = os.path.join(args.outdir, 'results.dump')
    if not os.path.isfile(saved):
        results, fasta, regions, stats = genResults(annot,kleats)
        writeFile(args.outdir, 'regions.bed', regions)
        writeFile(args.outdir, 'regions.fa', fasta)
        writeFile(args.outdir, 'stats', stats)
        pickle.dump(results, open(saved, 'wb'))
    else:
        results = pickle.load(open(saved, 'rb'))

    print 'DONE'

    ratios = computeRatios2(results, annot)

    ratios_path = os.path.join(args.outdir, 'ratios.dump')
    if not os.path.isfile(ratios_path):
        pickle.dump(ratios, open(ratios_path, 'wb'))
    else:
        ratios = pickle.load(open(ratios_path, 'rb'))

    ratios_human_path = os.path.join(args.outdir, 'ratios')
    with open(ratios_human_path, 'w') as f:
        for c in results:
            for g in results[c]:
                if g not in ratios:
                    continue
                for a in ratios[g]:
                    try:
                        f.write('{}\t{}\t{}\t{}\t{}\n'.format(c,g,annot[c][g][0].strand,a,(',').join([str(x) for x in ratios[g][a]])))
                    except KeyError:
                        print c, g
