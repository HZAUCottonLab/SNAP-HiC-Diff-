from check_intersect import check_intersect

def get_overlap(a, b):
        #print 'in get overlap'
        #print a
        #print b
        if a['chrom'] == b['chrom']:
                if (a['start']>b['start']) and (a['end']< b['end']): # community A is fully encompassed in B
                        #print "(a['start']>b['start']) and (a['end']< b['end'])"
                        overlap = a['end'] - a['start']
                        #print overlap
                        return overlap
                elif (b['start']>a['start']) and (b['end']< a['end']): # community B is fully encompassed in A
                        #print "(b['start']>a['start']) and (b['end']< a['end'])"
                        overlap = b['end'] - b['start']
                        #print overlap
                        return overlap
                elif check_intersect(a,b):

                        if a['start'] < b['end']:
                                #print "a['start'] < b['end']"
                                overlap = b['end'] - a['start']
                                #print overlap
                                return overlap
                        elif b['start'] < a['end']:
                                #print "b['start'] < a['end']"
                                overlap = a['end'] - b['start']
                                #print overlap
                        else:
                                print('intersect, but did not find an overlap!')
                else:
                        overlap = 0
                        #print overlap
                        return 0
        else:
                overlap = 0
                #print overlap
                return 0
