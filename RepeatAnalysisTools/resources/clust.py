import pysam,sys

NOCLUST = 999

def addHPtag(inBAM,outBAM,clusterMap,noCluster=NOCLUST,dropNoClust=False):
    '''clusterMap is map of {readname:cluster::int}'''
    with pysam.AlignmentFile(inBAM) as inbam:
        with pysam.AlignmentFile(outBAM,'wb',template=inbam) as outbam:
            for rec in inbam:
                if rec.query_name not in clusterMap and dropNoClust:
                    continue 
                clust = int(clusterMap.get(rec.query_name,noCluster))
                rec.set_tag('HP',clust)
                outbam.write(rec)
    return None

def readClusterFile(clustfile):
    res  = {}
    name = None
    cluster = 0
    with open(clustfile) as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:]
                cluster+=1
            else:
                #split off any extra fields
                read = '/'.join(line.split('/')[:3]) 
                res[read] = int(cluster)
    return res
