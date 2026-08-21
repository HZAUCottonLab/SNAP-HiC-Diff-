import numpy as np
from math import fsum
from scipy.stats import mode
# Updated 8/16/16 by HN such that the variance option is divided by N-1
def get_variance(Ci, bin_size):
        """
        Input: np array Ci of n_part partitions by communities (labelled 1 ... last community)
        Output: Ci after removing the community boundaries that are too great in variance
        V3 variances are computed off of a reiteration of nearby boundaries
        """
        n_part = np.size(Ci, 0) # number of partitions
        genomic_length = bin_size*np.size(Ci, 1) # length of partition in kb

        # Finds average of max community # along rows, then round it to find the average INTEGER community size
        avg_community_amount = int(round(np.mean(np.amax(Ci, 1))))

        # Empty 2D array to store the boundary indices between communities in each partition
        all_genomic_boundaries = []

        # Iterate to extract only the partitions with avg_community_amount communities
        for i in range(n_part):
                ar, current_indices = np.unique(Ci[i,:], return_index = True)
                current_indices.sort() # Sorted indices of boundaries

                # Change bins to genomic distance in kb
                genomic_boundaries = bin_size*current_indices

                # if number communities in current partition equals average community number
                if len(genomic_boundaries) == avg_community_amount:

                        # Concatenate the boundary indices except the first element (left boudndary)
                        all_genomic_boundaries.append(np.delete(genomic_boundaries, 0))

        # Most common boundary indices across all partitions after initial iteration
        if all([len(b)!=0 for b in all_genomic_boundaries]): # Only runs if all partitions have nonzero length
                temp_modes = mode(all_genomic_boundaries, 0)[0][0]

        else:
                print("ONLY 1 COMMUNITY EXISTS FOR THIS GAMMA")

        # Empty dictionary to store each column of boundary values
        boundary_dict = {}
        for i in range(avg_community_amount-1):
                boundary_dict[i] = []


        # Reiterate to include offset boundaries in variance calculation
        for i in range(n_part):
                failedcount = 0 # Number of missing boundaries at this partition

                ar, current_indices = np.unique(Ci[i,:], return_index = True)
                current_indices.sort() # Sorted indices of boundaries

                # Remove the first boundary, 0
                current_indices = np.delete(current_indices, 0)

                # Array of genomic coordinates in kb
                genomic_boundaries = bin_size*current_indices

                # Append any additional boundaries to consider
                for j in range(avg_community_amount-1): # Iterate over every boundary mode

                        offsets = genomic_boundaries - temp_modes[j] # Difference each boundary is to the most common boundary
                        within_50 = [] # Temporary list that contains all the boundaries within 50% community to the current temp_mode[j]

                        # size of community to left of boundary
                        if j == 0:
                                community_left = temp_modes[j] - 0
                        else:
                                community_left = temp_modes[j] - temp_modes[j-1]

                        # size of community to right of boundary
                        if j == len(temp_modes) - 1:
                                community_right = genomic_length - temp_modes[j]
                        else:
                                community_right = temp_modes[j+1] - temp_modes[j]

                        # Checks if boundary to be compared is on the left side, right side, or right on the mode
                        for key, offset in enumerate(offsets):
                                if (offset < 0 and abs(offset) < 0.5*(community_left)) or \
                                (offset > 0 and abs(offset) < 0.5*(community_right)) or offset == 0:
                                        within_50.append(genomic_boundaries[key]) # Append the genomic boundary that is within 50% community of boundary

                        # If not empty
                        if within_50:
                                # Out of within_50, the closest one to the mode boundary
                                closest = within_50[abs(within_50 - temp_modes[j]).argmin()]

                                # Include it into the boundary dictionary to be considered
                                boundary_dict[j].append(closest)
                        else:
                                failedcount+=1

                #print "%d out of %d modes at this partition did not find a real boundary" %(failedcount, avg_community_amount-1)
                #print "This partition has %d boundaries\n" %len(genomic_boundaries)

        # Means and variances of each boundary after including the offset partitions
        boundary_means = []
        boundary_vars = []
        for i in range(avg_community_amount-1):
                #print "boundary across partitions for %d:\n" %i
                #print boundary_dict[i]
                boundary_means.append(np.mean(boundary_dict[i]))
                boundary_vars.append(np.var(boundary_dict[i]))

        # Inserts and appends the first and last values
        boundary_means = np.array(boundary_means)
        boundary_vars = np.array(boundary_vars)

        boundary_means = np.append(boundary_means, genomic_length)
        boundary_vars = np.append(boundary_vars, 0)
        boundary_means = np.insert(boundary_means, 0, 0)
        boundary_vars = np.insert(boundary_vars, 0, 0)
        #print "boundary_means: ", boundary_means
        #print "boundary_vars: ", boundary_vars

        # Return non-np lists so that the text does not wrap
        return (list(boundary_means), list(boundary_vars))


def get_variance_about_boundaries_size_thresholded(Ci, markers, size_threshold, return_edges=True, type='regular'):
        """ Determines the variance within a partition block about input marker boundaries

        Inputs
        ------
                Ci: ndarray
                        partition block
                markers: 1xn list or array
                        boundary markers indexing Ci about which real boundaries should be searched for, should not contain 0 or the length of Ci
                size_threshold: integer
                        number of nodes a community must be larger than to be considered in the variance computation
                return_edges: boolean
                        if true, appends 0 to beginning and end of boundary_vars to match the left and right walls of Ci
                type: str
                        'regular' : variances are computed as sum((boundaries - marker)^2)/N
                        'mean' : variances are computed as sum(abs(boundaries-marker))/N, the mean distance from marker
                        'mode' : variances are computed as mode(abs(boundaries-marker)), the most common distance from markers

        Outputs
        -------
                boundary_vars: list of floats
                        variances corresponding to each element in markers.
                distances: nested list of integers
                        distances from each boundary to its closest marker. Does not contain elements for left and right edges of partition block
        """


        # number of partitions
        n_part = np.size(Ci, 0)
        num_nodes = np.size(Ci, 1)
        markers.sort()

        # number of nodes
        node_length = np.size(Ci, 1)

        # key = marker, value = list of boundaries found around marker
        boundary_dict = {}
        previous_marker = 0
        markers_size_thresholded = []
        for i in range(len(markers)):
                boundary_dict[markers[i]] = []
        #for marker in markers:
                if i < len(markers) - 1:
                        if abs(markers[i] - previous_marker) > size_threshold or abs(markers[i] - markers[i+1]) > size_threshold:
                                markers_size_thresholded.append(markers[i])
                elif i == len(markers) - 1:
                        if abs(markers[i] - previous_marker) > size_threshold or abs(markers[i] - num_nodes) > size_threshold:
                                markers_size_thresholded.append(markers[i])
                previous_marker = markers[i]

        #print 'length of markers', len(markers)
        #print 'length of size-thresholded markers', len(markers_size_thresholded)
        #print 'number of nodes', num_nodes
        #print 'markers', markers
        #print 'markers-size_thresholded', markers_size_thresholded

        for i in range(n_part):
                # Number of times a given partition fails to find a boundary proximal to a marker
                failedcount = 0

                # Indices within the current partition that mark the beginning of each community
                current_boundaries = np.where(Ci[i][1:] != Ci[i][:-1])[0]+1

                # Iterate over each marker boundary
                for j in range(len(markers)):

                        # displacement of each partition boundary to the marker
                        offsets = current_boundaries - markers[j]

                        # Temporary list that contains all the boundaries within 50% community size to the current marker
                        within_50 = []

                        # size of community to left of marker
                        if j == 0:
                                community_left = markers[j] - 0
                        else:
                                community_left = markers[j] - markers[j-1]

                        # size of community to right of boundary
                        if j == len(markers) - 1:
                                community_right = node_length - markers[j]
                        else:
                                community_right = markers[j+1] - markers[j]

                        # Add boundaries that are within 50% of community size to the current marker
                        for key, offset in enumerate(offsets):
                                if (offset < 0 and abs(offset) < 0.5*(community_left)) or \
                                (offset > 0 and abs(offset) < 0.5*(community_right)) or offset == 0:
                                        within_50.append(current_boundaries[key])

                        # If any boundaries are appended to within_50
                        if within_50:

                                # Out of within_50, the closest boundary to the marker
                                closest = within_50[abs(within_50 - markers[j]).argmin()]

                                # Include it into the boundary dictionary to be considered at the current marker
                                boundary_dict[markers[j]].append(closest)

                        else:

                                # This partition did not find any boundaries around the current marker
                                failedcount+=1


        # each element is a nested list of distances from boundaries close to a marker to that marker
        distances = []
        distances_squared = []
        # each element is a variance of the boundaries close to a marker about that expected marker
        boundary_vars = []

        for marker in markers:
                if marker in markers_size_thresholded:
                        #print 'marker in size_thresholded', marker
                #print 'how many boundaries around marker %d: ' %marker, len(boundary_dict[marker])
                        distances.append([abs(bound - marker) for bound in boundary_dict[marker]])
                        distances_squared.append([(bound-marker)**2 for bound in boundary_dict[marker]])
                        boundary_vars.append(var_about_center(boundary_dict[marker], marker, type=type))

        # add 0s to both side of boundary_vars corresponding to the left and right walls of Ci
        if return_edges:
                boundary_vars.append(0)
                boundary_vars.insert(0, 0)

        return boundary_vars, distances, distances_squared, markers_size_thresholded



def get_variance_about_boundaries(Ci, markers, return_edges=True, type='regular'):
        """ Determines the variance within a partition block about input marker boundaries

        Inputs
        ------
                Ci: ndarray
                        partition block
                markers: 1xn list or array
                        boundary markers indexing Ci about which real boundaries should be searched for, should not contain 0 or the length of Ci
                return_edges: boolean
                        if true, appends 0 to beginning and end of boundary_vars to match the left and right walls of Ci
                type: str
                        'regular' : variances are computed as sum((boundaries - marker)^2)/N
                        'mean' : variances are computed as sum(abs(boundaries-marker))/N, the mean distance from marker
                        'mode' : variances are computed as mode(abs(boundaries-marker)), the most common distance from markers

        Outputs
        -------
                boundary_vars: list of floats
                        variances corresponding to each element in markers.
                distances: nested list of integers
                        distances from each boundary to its closest marker. Does not contain elements for left and right edges of partition block
        """

        # number of partitions
        n_part = np.size(Ci, 0)

        # number of nodes
        node_length = np.size(Ci, 1)

        # key = marker, value = list of boundaries found around marker
        boundary_dict = {}
        for marker in markers:
                boundary_dict[marker] = []

        for i in range(n_part):
                # Number of times a given partition fails to find a boundary proximal to a marker
                failedcount = 0

                # Indices within the current partition that mark the beginning of each community
                current_boundaries = np.where(Ci[i][1:] != Ci[i][:-1])[0]+1

                # Iterate over each marker boundary
                for j in range(len(markers)):

                        # displacement of each partition boundary to the marker
                        offsets = current_boundaries - markers[j]

                        # Temporary list that contains all the boundaries within 50% community size to the current marker
                        within_50 = []

                        # size of community to left of marker
                        if j == 0:
                                community_left = markers[j] - 0
                        else:
                                community_left = markers[j] - markers[j-1]

                        # size of community to right of boundary
                        if j == len(markers) - 1:
                                community_right = node_length - markers[j]
                        else:
                                community_right = markers[j+1] - markers[j]

                        # Add boundaries that are within 50% of community size to the current marker
                        for key, offset in enumerate(offsets):
                                if (offset < 0 and abs(offset) < 0.5*(community_left)) or \
                                (offset > 0 and abs(offset) < 0.5*(community_right)) or offset == 0:
                                        within_50.append(current_boundaries[key])

                        # If any boundaries are appended to within_50
                        if within_50:

                                # Out of within_50, the closest boundary to the marker
                                closest = within_50[abs(within_50 - markers[j]).argmin()]

                                # Include it into the boundary dictionary to be considered at the current marker
                                boundary_dict[markers[j]].append(closest)

                        else:

                                # This partition did not find any boundaries around the current marker
                                failedcount+=1

                #print "%d out of %d markers failed to find a boundary within this partition" %(failedcount, len(markers))
                #print "This partition has %d boundaries\n" %len(current_boundaries)

        # each element is a nested list of distances from boundaries close to a marker to that marker
        distances = []
        distances_squared = []
        # each element is a variance of the boundaries close to a marker about that expected marker
        boundary_vars = []

        for marker in markers:
                #print 'how many boundaries around marker %d: ' %marker, len(boundary_dict[marker])
                distances.append([abs(bound - marker) for bound in boundary_dict[marker]])
                distances_squared.append([(bound-marker)**2 for bound in boundary_dict[marker]])
                boundary_vars.append(var_about_center(boundary_dict[marker], marker, type=type))

        # add 0s to both side of boundary_vars corresponding to the left and right walls of Ci
        if return_edges:
                boundary_vars.append(0)
                boundary_vars.insert(0, 0)

        return boundary_vars, distances, distances_squared

def var_about_center(nums, center, type='regular'):
        """ takes the variance of a set of numbers about a specified center number

        Inputs
        ------
                nums: list of int/float
                        numbers to calculate variance on
                center: int/float
                        expected value to find variance around
                type: str
                        type of variance to return.
                        'regular' = sum((X - X-bar)^2)/(N-1)
                        'mean' = sum(abs(X - X-bar))/N
                        'mode' = mode(abs(X - X-bar))

        Outputs
        -------
                var: float
                        variance of nums
        """

        # Return nan if nums is empty
        if len(nums) <2:
                print("NO VALUES TO FIND VARIANCE FOR")
                return float('nan')
        #if not nums:
                #return float('nan')

        # warning if center is outside of range of nums
        if center > max(nums) or center < min(nums):
                print("warning: center is outside of values")

        # type of variance to return
        if type == 'regular':
                var = fsum([(float(x) - float(center))**2 for x in nums])/(float(len(nums))-1)
        elif type == 'mean':
                var = fsum([abs(float(x) - float(center)) for x in nums])/float(len(nums))
        elif type == 'mode':
                var = mode([abs(float(x) - float(center)) for x in nums])[0][0]

        return var
