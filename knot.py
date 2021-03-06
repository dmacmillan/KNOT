
__version__ = '1.0.2'

import argparse, sys, os, math
import colorsys
import cPickle as pickle
from copy import deepcopy
import pysam
import time
from Queue import Queue
from threading import Thread
import subprocess

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

def clsGTF(grouped):
    classes = {'i': set(),
               'ii': set(),
               'iii': set()}
    for c in grouped:
        intervals = []
        for g in grouped[c]:
            #for gtf in grouped[c][g]:
            strand = grouped[c][g][0].strand
            if strand == '-':
                grouped[c][g] = grouped[c][g][::-1]
            intervals.append(grouped[c][g][-1])
        intervals.sort(key=lambda x: x.start)
        for i in xrange(len(intervals)-1):
            for j in xrange(i+1, len(intervals)):
                if intervals[j].start > intervals[i].end:
                    continue
                classes['ii'].add(intervals[j].attribute['gene_name'])
                classes['ii'].add(intervals[i].attribute['gene_name'])
    return classes

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

def kleatLinkage(sites, window=20):#atuple):
    #sites, window = atuple
    length = len(sites)
    #print length
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
            kleatLinkage(sites, window)#(sites, window))
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
        result.append('color="{}"'.format((',').join([str(x) for x in color])))
    return 'track ' + (' ').join(result)

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

def computeRatios(results, annot):
    ratios = {}
    output = {}
    for chrom in results:
        output[chrom] = {}
        for gene in results[chrom]:
            output[chrom][gene] = {}
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
                        ratios[gene][align] = {'med': [results[chrom][gene][span][align]['med']],'avg': [results[chrom][gene][span][align]['avg']]}
                    else:
                        ratios[gene][align]['med'].append(results[chrom][gene][span][align]['med'])
                        ratios[gene][align]['avg'].append(results[chrom][gene][span][align]['avg'])
            samples = output[chrom][gene]
            for s in ratios[gene]:
                length = len(ratios[gene][s]['med'])
                if length <= 1:
                    continue
                for i in xrange(length-1):
                    for j in xrange(i+1,length):
                        med_i = ratios[gene][s]['med'][i]
                        avg_i = ratios[gene][s]['avg'][i]
                        med_j = ratios[gene][s]['med'][j]
                        avg_j = ratios[gene][s]['avg'][j]
                        avg_ratio = float(avg_j)/avg_i
                        med_ratio = float(med_j)/med_i
                        if s not in samples:
                            samples[s] = {'avg_ratios':[avg_ratio],
                                          'med_ratios':[med_ratio],
                                          'regions': keys}
                        else:
                            samples[s]['avg_ratios'].append(avg_ratio)
                            samples[s]['med_ratios'].append(med_ratio)
            #output.append(samples)
    return output

def interpretRatios(ratios):
    diffs = []
    for chrom in ratios:
        for gene in ratios[chrom]:
            c = ratios[chrom][gene]
            if not c:
                continue
            if len(set([len(c[s]['med_ratios']) for s in c])) != 1:
                continue
            samples = c.keys()
            for i in xrange(len(samples)-1):
                for j in xrange(i+1,len(samples)):
                    for reg in xrange(len(c[samples[i]]['regions'])-1):
                        try:
                            avg_ratios_i = c[samples[i]]['avg_ratios'][reg]
                            avg_ratios_j = c[samples[j]]['avg_ratios'][reg]
                            med_ratios_i = c[samples[i]]['med_ratios'][reg]
                            med_ratios_j = c[samples[j]]['med_ratios'][reg]
                            avg_diff = abs(avg_ratios_i - avg_ratios_j)
                            med_diff = abs(med_ratios_i - med_ratios_j)
                            thing = {'chrom': chrom,
                                     'gene': gene,
                                     'sample_i': samples[i],
                                     'sample_j': samples[j],
                                     'avg_diff': avg_diff,
                                     'med_diff': med_diff,
                                     'region': c[samples[i]]['regions'][reg]}
                            diffs.append(thing)
                        except IndexError:
                            continue
    return diffs
                        

def interpretRatios(data):
    diffs = []
    for item in data:
        samples = item.keys()
        for i in xrange(len(samples)-1):
            for j in xrange(i+1,len(samples)):
                for reg in xrange(len(item[samples[i]]['regions'])-1):
                    try:
                        avg_ratios_i = item[samples[i]]['avg_ratios'][reg]
                        avg_ratios_j = item[samples[j]]['avg_ratios'][reg]
                        med_ratios_i = item[samples[i]]['med_ratios'][reg]
                        med_ratios_j = item[samples[j]]['med_ratios'][reg]
                        avg_diff = abs(avg_ratios_i - avg_ratios_j)
                        med_diff = abs(med_ratios_i - med_ratios_j)
                        thing = {'chrom': item[samples[i]]['chrom'],
                                 'gene': item[samples[i]]['gene'],
                                 'sample_i': samples[i],
                                 'sample_j': samples[j],
                                 'avg_diff': avg_diff,
                                 'med_diff': med_diff,
                                 'region': item[samples[i]]['regions'][reg]}
                        diffs.append(thing)
                    except IndexError:
                        continue
    #avg_stdev = np.std([x['avg_diff'] for x in diffs],ddof=1)
    #med_stdev = np.std([x['med_diff'] for x in diffs],ddof=1)
    return diffs

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

def genRegions(annot, kleats, cls):
    results = {}
    regions = ''
    for chrom in annot:
        if chrom not in kleats:
            continue
        results[chrom] = {}
        for gene in annot[chrom]:
            if gene not in kleats[chrom]:
                continue
            if (not annot[chrom][gene]):
                continue
            if gene in cls['ii']:
                continue
            results[chrom][gene] = {}
            strand = annot[chrom][gene][0].strand
            if strand == '-':
                keys = annot[chrom][gene][:1]
            else:
                keys = annot[chrom][gene][-2:]
            for region in keys:
                last = region.start
                cleaved = False
                cuts = []
                for k in kleats[chrom][gene]:
                    if (region.start < k < region.end):
                        cleaved = True
                        cuts.append(k)
                if not cleaved:
                    continue
                coords = sorted([region.start, region.end] + cuts)
                my_regions = []
                for i in xrange(len(coords)-1):
                    start = coords[i]
                    stop = coords[i+1]
                    if abs(stop-start) <= 20:
                        coords[i+1] = coords[i]
                        continue
                    my_regions.append([start,stop])
                for i,r in enumerate(my_regions):
                    regions += ('\t').join([chrom, str(r[0]), str(r[1]), '{}_utr_{}'.format(gene,i), '0', region.strand, '\n'])
                    span = '{}-{}'.format(r[0],r[1])
                    results[chrom][gene][span] = {}
    return results, regions

def getPileupRegion(results, aligns, sample, chrom, gene, start, stop):
    results, aligns, sample, chrom, gene, start, stop
    print 'Getting pileup for {} {} {}:{}-{}'.format(sample, gene, chrom, start, stop)
    r = []
    bg = aligns[sample]['bg']
    for p in aligns[sample]['align'].pileup(chrom, start, stop):
        r.append(p.n)
        bg.append('{}\t{}\t{}\t{}'.format(chrom,p.pos,p.pos+1,p.n))
    length = len(r)
    if length < 1:
        return
    r = sorted(r)
    median = calcMedian(r)
    mean = sum(r)/length
    _min = r[0]
    _max = r[-1]
    q1, q3 = calcIQR(r)
    se = sstdev(r)/math.sqrt(length)
    span = '{}-{}'.format(start,stop)
    results[chrom][gene][span][sample] = {'median': median, 'mean': mean, 'q1': q1, 'q3': q3, 'min': _min, 'max': _max, 'se': se}
    return

def genStats(results, aligns, annot):
    data = [('\t').join(['GENE','STRAND','SAMPLE','REGION','LENGTH','MIN','Q1','MED','Q3','MAX','MEAN','SE'])]
    for chrom in results:
        for gene in results[chrom]:
            c = results[chrom][gene]
            if not c:
                continue
            strand = annot[chrom][gene][0].strand
            for span in c:
                start, stop = [int(x) for x in span.split('-')]
                for sample in c[span]:
                    cur = c[span][sample]
                    data.append(('\t').join([str(x) for x in [gene, strand, sample, chrom + ':' + span, stop-start, cur['min'], cur['q1'], cur['median'], cur['q3'], cur['max'], cur['mean'], cur['se']]]))
    return data

def doEverything(annot, kleats, cls, myfile):
    data = ('\t').join(['GENE','STRAND','SAMPLE','REGION','LENGTH','MIN','Q1','MED','Q3','MAX','MEAN','SE']) + '\n'
    myfile.write(data)
    myfile.flush()
    for chrom in kleats:
        #print 'chrom: {}'.format(chrom)
        for gene in kleats[chrom]:
            #print '  gene: {}'.format(gene)
            if not kleats[chrom][gene]:
                #print 'kleats {} {} is empty'.format(chrom,gene)
                continue
            elif gene in cls['ii']:
                #print '{} is overlapping'.format(gene)
                continue
            elif chrom not in annot:
                #print '{} is not in annot'.format(chrom)
                continue
            elif gene not in annot[chrom]:
                #print '{} not in annot {}'.format(gene, chrom)
                continue
            strand = annot[chrom][gene][0].strand
            if strand == '-':
                utr3 = annot[chrom][gene][:1]
            else:
                utr3 = annot[chrom][gene][-2:]
            for region in utr3:
                #print '    span: {}-{}'.format(region.start, region.end)
                last = region.start
                intervals = []
                splices = 0
                cleaved = False
                cuts = []
                for k in kleats[chrom][gene]:
                    if (region.start < k < region.end):
                        cleaved = True
                        cuts.append(k)
                if not cleaved:
                    continue
                coords = sorted([region.start, region.end] + cuts)
                my_regions = []
                for i in xrange(len(coords)-1):
                    start = coords[i]
                    stop = coords[i+1]
                    if abs(stop-start) <= 20:
                        coords[i+1] = coords[i]
                        continue
                    my_regions.append([start,stop])
                for i,r in enumerate(my_regions):
                    span = '{}-{}'.format(r[0],r[1])
                    for sample in aligns:
                        s = []
                        #bg = aligns[sample]['bg']
                        for p in aligns[sample]['align'].pileup(chrom, r[0], r[1]):
                            if r[0] <= p.pos <= r[1]:
                                s.append(p.n)
                            #bg.append('{}\t{}\t{}\t{}'.format(chrom,p.pos,p.pos+1,p.n))
                        length = len(s)
                        if length < 3:
                            continue
                        s = sorted(s)
                        median = calcMedian(s)
                        mean = sum(s)/length
                        _min = s[0]
                        _max = s[-1]
                        q1, q3 = calcIQR(s)
                        se = sstdev(s)/math.sqrt(length)
                        span = '{}-{}'.format(r[0],r[1])
                        print 'Writing {}\t{}\t{}'.format(sample, gene, chrom+':'+span)
                        towrite = ('\t').join([str(x) for x in [gene, strand, sample, chrom + ':' + span, r[1]-r[0], _min, q1, median, q3, _max, mean, se]]) + '\n'
                        myfile.write(towrite)
                        myfile.flush()

def genResults(results, annot, kleats, cls):
    q = Queue()
    for chrom in annot:
        if chrom not in kleats:
            continue
        for gene in annot[chrom]:
            if gene not in kleats[chrom]:
                continue
            if (not annot[chrom][gene]):
                continue
            if gene in cls['ii']:
                continue
            strand = annot[chrom][gene][0].strand
            if strand == '-':
                annot[chrom][gene] = annot[chrom][gene][:1]
            else:
                annot[chrom][gene] = annot[chrom][gene][-2:]
            for region in annot[chrom][gene]:
                last = region.start
                intervals = []
                splices = 0
                cleaved = False
                cuts = []
                for k in kleats[chrom][gene]:
                    if (region.start < k < region.end):
                        cleaved = True
                        cuts.append(k)
                if not cleaved:
                    continue
                coords = sorted([region.start, region.end] + cuts)
                my_regions = []
                for i in xrange(len(coords)-1):
                    start = coords[i]
                    stop = coords[i+1]
                    if abs(stop-start) <= 20:
                        coords[i+1] = coords[i]
                        continue
                    my_regions.append([start,stop])
                temp = []
                for i,r in enumerate(my_regions):
                    span = '{}-{}'.format(r[0],r[1])
                    results[chrom][gene][span] = {}
                    for sample in aligns:
                        getPileupRegion(results, aligns, sample, chrom, gene, r[0], r[1])
                        #q.put([results, aligns, sample, chrom, gene, r[0], r[1]])
    return

def calcMedian(lst):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    if lstLen == 0:
        return None
    elif lstLen == 1:
        return lst[0]
    index = (lstLen - 1) // 2
    if (lstLen % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1])/2.0

def sprint(text):
    sys.stdout.write(text)
    sys.stdout.flush()

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

def worker(q):
    while True:
        item = q.get()
        getPileupRegion(*item)
        q.task_done()

def AHC(mylist, window=20):
    length = len(mylist)
    if length <= 1:
        return mylist
    mymin = float('inf')
    r = s = None
    for i in xrange(length-1):
        dist = mylist[i+1] - mylist[i]
        if dist < mymin:
            r = i
            s = i + 1
            mymin = dist
    if mymin <= window:
        mylist[r] = sum([mylist[r],mylist[s]])/2
        del(mylist[s])
        AHC(mylist, window)
    return mylist

def checkSeq(align, region, samtools):
    command = [samtools, 'view', align, region]
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    o,e = p.communicate()
    if e == 'Container header CRC32 failure\n':
        return False
    return True

# Begin main thread
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Compare two bedgraph tracks")
    
    parser.add_argument('annotation', help='Genome annotation file if GTF format')
    parser.add_argument('kleats', help='File containing list of kleat output files to use')
    parser.add_argument('alignments', help='File containing list of alignment output files to use. Must be in BAM or SAM format. Names of each BAM/SAM file must also differ')
    parser.add_argument('-a', '--align', nargs='*', help='One or more alignment files to process, overrides the "alignments parameter"')
    parser.add_argument('-cw', '--cluster_window', type=int, default=20, help='Set the window size for clustering KLEAT cleavage sites. Default = 20')
    parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Directory to output to. Default is current directory')
    parser.add_argument('-r', '--reference', default='/home/dmacmillan/references/hg19/hg19.fa', help='Path to the reference genome from which to fetch sequences')
    parser.add_argument('-la', '--load_annotation', help='annotation.dump file')
    parser.add_argument('-lc', '--load_clusters', help='clusters.dump file')
    
    args = parser.parse_args()

    if not os.path.isdir(args.outdir):
        print 'Outdir {} does not exist, creating ...'.format(args.outdir)
        os.makedirs(args.outdir)

    # Output files
    results_dump = os.path.join(args.outdir, 'results.dump')
    annot_dump = os.path.join(args.outdir, 'annotation.dump')
    clusters_dump = os.path.join(args.outdir, 'clusters.dump')

    ref = pysam.FastaFile(args.reference)

    samtools = '/gsc/btl/linuxbrew/bin/samtools'

    all_kleats = []

    aligns = {}

    if not args.align:
        num_aligns = 0
        with open(args.alignments, 'r') as f:
            for line in f:
                num_aligns += 1
        for i,b in enumerate(parseConfig(args.alignments)):
            sprint('Loading alignment {}/{}\r'.format(i,num_aligns))
            name,ftype = os.path.splitext(os.path.abspath(b))
            name = os.path.basename(name)
            if ftype == '.cram':
                aligns[name] = {'align': pysam.AlignmentFile(b, 'rc'),'bg': None, 'path': b}
            else:
                aligns[name] = {'align': pysam.AlignmentFile(b, 'rb'),'bg': None, 'path': b}
        print 'Loading alignment {}/{} DONE'.format(num_aligns,num_aligns)
    else:
        for i,b in enumerate(args.align):
            name,ftype = os.path.splitext(os.path.abspath(b))
            name = os.path.basename(name)
            if ftype == '.cram':
                aligns[name] = {'align': pysam.AlignmentFile(b, 'rc'),'bg': None, 'path': b}
            else:
                aligns[name] = {'align': pysam.AlignmentFile(b, 'rb'),'bg': None, 'path': b}
            if not aligns[name]['align'].has_index():
                sprint('{} has no index! Creating... '.format(name))
                pysam.index(b)
                aligns[name] = {'align': pysam.AlignmentFile(b, 'rc'),'bg': None, 'path': b}
                print 'DONE'

    N = len(aligns)
    HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    RGB_tuples = map(lambda x: [255*y for y in x],RGB_tuples)

    for i,name in enumerate(aligns):
        aligns[name]['bg'] = [genTrackLine(name+'.bg', name+'.bg', 'bedGraph', color=RGB_tuples[i])]

    bad_aligns = []
    quit = False
    # Check alignments
    for sample in aligns:
        if not checkSeq(aligns[sample]['path'], 'chr1:10000-15000', samtools):
            quit = True
            print 'BAD ALIGNMENT: {}'.format(sample)
            bad_aligns.append('/gsc/btl/linuxbrew/bin/samtools view -@ 12 -b -h /projects/btl/polya/sorted_star_ccle/{}.cram | /gsc/btl/linuxbrew/bin/samtools view -C -T /home/dmacmillan/references/hg19/hg19.fa -@ 12 -o /projects/btl/polya/sorted_star_ccle/reconverted_{}.cram -; mv /projects/btl/polya/sorted_star_ccle/{}.cram /projects/btl/polya/sorted_star_ccle/old_{}.cram; rm /projects/btl/polya/sorted_star_ccle/{}.cram.crai; mv /projects/btl/polya/sorted_star_ccle/reconverted_{}.cram /projects/btl/polya/sorted_star_ccle/{}.cram; /gsc/btl/linuxbrew/bin/samtools index /projects/btl/polya/sorted_star_ccle/{}.cram &'.format(sample, sample,  sample, sample, sample, sample, sample, sample))

    if quit:
        writeFile(args.outdir, 'fix_bad_aligns', *bad_aligns)
        sys.exit('Some alignments were broken, please fix before continuing')

    if not args.load_clusters:
        # Parse KLEAT files
        for i,k in enumerate(parseConfig(args.kleats)):
            sprint('Loading kleat {}/{}\r'.format(i,num_aligns))
            kleat = parseKleat(k, with_pas=True, min_num_bridge_reads=2, min_bridge_read_tail_len=4)
            all_kleats += kleat
        print 'Loading kleat {}/{} DONE'.format(num_aligns,num_aligns)
        
        # Grouping
        sprint('Grouping kleat data ...')
        kleats = groupKleat(all_kleats)
        print 'DONE'
    
        # Clustering
        start = time.time()
        sprint('Clustering cleavage events ...')
        for chrom in kleats:
            for gene in kleats[chrom]:
                if not kleats[chrom][gene]:
                    continue
                kleats[chrom][gene] = sorted(set([x.cleavage_site for x in kleats[chrom][gene]]))
                kleats[chrom][gene] = AHC(kleats[chrom][gene], args.cluster_window)
        print 'DONE ({}s)'.format(time.time()-start)
    
        pickle.dump(kleats, open(clusters_dump, 'wb'))
    else:
        sprint('Loading clusters ...')
        kleats = pickle.load(open(args.load_clusters, 'rb'))
        print 'DONE'

    if not args.load_annotation:
        sprint('Parsing GTF ...')
        annot = parseGTF(args.annotation, seqnames=['chr{}'.format(x) for x in range(1,23)] + ['chrX', 'chrY'], sources=['protein_coding'], features='UTR')
    
        annot = groupGTF(annot)
        pickle.dump(annot, open(annot_dump, 'wb'))
        cls = clsGTF(annot)
        print 'DONE'
    else:
        sprint('Loading annotation ...')
        annot = pickle.load(open(args.load_annotation,'rb'))
        cls = clsGTF(annot)
        print 'DONE'

    myfile = open(os.path.join(args.outdir, 'stats'), 'w')

    doEverything(annot, kleats, cls, myfile)

    myfile.close()
    #results, regions = genRegions(annot, kleats, cls)

    #writeFile(args.outdir, 'regions.bed', regions)

    #q = genResults(results, annot, kleats, cls)
    #print q.qsize()

#    for i in range(args.threads):
#        t = Thread(target = worker, args = (q, ))
#        t.setDaemon(True)
#        t.start()
#
#    q.join()
    #print 'DONE'

    #sprint('Writing bedgraphs ...')
    #for a in aligns:
    #    writeFile(args.outdir, a+'.bg', *aligns[a]['bg'])
    #print 'DONE'

    #sprint('Generating statistics ...')
    #stats = genStats(results, aligns, annot)
    #print 'DONE'

    #writeFile(args.outdir, 'stats', *stats)

    #ratios = computeRatios(results, annot)
    #diffs = interpretRatios(ratios)
