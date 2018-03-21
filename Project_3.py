"""implementation of project 3 for algorithmic thinking"""
import math


class Cluster:
    """
    Class for creating and merging clusters of counties
    """

    def __init__(self, fips_codes, horiz_pos, vert_pos, population, risk):
        """
        Create a cluster based the models a set of counties' data
        """
        self._fips_codes = fips_codes
        self._horiz_center = horiz_pos
        self._vert_center = vert_pos
        self._total_population = population
        self._averaged_risk = risk

    def __repr__(self):
        """
        String representation assuming the module is "alg_cluster".
        """
        rep = "Cluster("
        rep += str(self._fips_codes) + ", "
        rep += str(self._horiz_center) + ", "
        rep += str(self._vert_center) + ", "
        rep += str(self._total_population) + ", "
        rep += str(self._averaged_risk) + ")"
        return rep

    def fips_codes(self):
        """
        Get the cluster's set of FIPS codes
        """
        return self._fips_codes

    def horiz_center(self):
        """
        Get the averged horizontal center of cluster
        """
        return self._horiz_center

    def vert_center(self):
        """
        Get the averaged vertical center of the cluster
        """
        return self._vert_center

    def total_population(self):
        """
        Get the total population for the cluster
        """
        return self._total_population

    def averaged_risk(self):
        """
        Get the averaged risk for the cluster
        """
        return self._averaged_risk

    def copy(self):
        """
        Return a copy of a cluster
        """
        copy_cluster = Cluster(set(self._fips_codes), self._horiz_center, self._vert_center,
                               self._total_population, self._averaged_risk)
        return copy_cluster

    def distance(self, other_cluster):
        """
        Compute the Euclidean distance between two clusters
        """
        vert_dist = self._vert_center - other_cluster.vert_center()
        horiz_dist = self._horiz_center - other_cluster.horiz_center()
        return math.sqrt(vert_dist ** 2 + horiz_dist ** 2)

    def merge_clusters(self, other_cluster):
        """
        Merge one cluster into another
        The merge uses the relatively populations of each
        cluster in computing a new center and risk

        Note that this method mutates self
        """
        if len(other_cluster.fips_codes()) == 0:
            return self
        else:
            self._fips_codes.update(set(other_cluster.fips_codes()))

            # compute weights for averaging
            self_weight = float(self._total_population)
            other_weight = float(other_cluster.total_population())
            self._total_population = self._total_population + other_cluster.total_population()
            self_weight /= self._total_population
            other_weight /= self._total_population

            # update center and risk using weights
            self._vert_center = self_weight * self._vert_center + other_weight * other_cluster.vert_center()
            self._horiz_center = self_weight * self._horiz_center + other_weight * other_cluster.horiz_center()
            self._averaged_risk = self_weight * self._averaged_risk + other_weight * other_cluster.averaged_risk()
            return self

    def cluster_error(self, data_table):
        """
        Input: data_table is the original table of cancer data used in creating the cluster.

        Output: The error as the sum of the square of the distance from each county
        in the cluster to the cluster center (weighted by its population)
        """
        # Build hash table to accelerate error computation
        fips_to_line = {}
        for line_idx in range(len(data_table)):
            line = data_table[line_idx]
            fips_to_line[line[0]] = line_idx

        # compute error as weighted squared distance from counties to cluster center
        total_error = 0
        counties = self.fips_codes()
        for county in counties:
            line = data_table[fips_to_line[county]]
            singleton_cluster = Cluster(set([line[0]]), line[1], line[2], line[3], line[4])
            singleton_distance = self.distance(singleton_cluster)
            total_error += (singleton_distance ** 2) * singleton_cluster.total_population()
        return total_error


def slow_closest_pair(cluster_list):
    """Slow way of finding closest pair of clusters"""
    (dist, idx1, idx2) = (float('Inf'), -1, -1)
    for cluster1 in cluster_list:
        for cluster2 in cluster_list:
            if cluster_list.index(cluster1) != cluster_list.index(cluster2):
                (dist, idx1, idx2) = min(tuple([dist, idx1, idx2]),
                                         tuple([cluster1.distance(cluster2), cluster_list.index(cluster1),
                                                cluster_list.index(cluster2)]), key=lambda tup: tup[0])
    return tuple([dist, idx1, idx2])


def closest_pair_strip(cluster_list, horiz_center, half_width):
    """helper function, no idea what it really does"""
    clusters = [index for index in range(len(cluster_list))
                if abs(cluster_list[index].horiz_center() - horiz_center) < half_width]
    clusters.sort(key=lambda idx: cluster_list[idx].vert_center())
    num = len(clusters)
    (dist, idx1, idx2) = (float('Inf'), -1, -1)
    for thing in range(0, num - 1):
        for thing2 in range(thing + 1, min(thing + 3, num - 1) + 1):
            (dist, idx1, idx2) = min((dist, idx1, idx2),
                                     (cluster_list[clusters[thing]].distance(cluster_list[clusters[thing2]]),
                                      clusters[thing], clusters[thing2]), key=lambda tup: tup[0])
    if idx2 < idx1:
        return tuple([dist, idx2, idx1])
    return tuple([dist, idx1, idx2])


def fast_closest_pair(cluster_list):
    """fast way to find closest pair of clusters"""
    num = len(cluster_list)
    if num <= 3:
        (dist, idx1, idx2) = slow_closest_pair(cluster_list)
    else:
        lst_m = num / 2
        cluster_left = cluster_list[:lst_m]
        cluster_right = cluster_list[lst_m:]
        (distl, idx1l, idx2l) = fast_closest_pair(cluster_left)
        (distr, idx1r, idx2r) = fast_closest_pair(cluster_right)
        (dist, idx1, idx2) = min((distl, idx1l, idx2l), (distr, idx1r + lst_m, idx2r + lst_m), key=lambda tup: tup[0])
        mid = 0.5 * (cluster_list[lst_m - 1].horiz_center() + cluster_list[lst_m].horiz_center())
        (dist, idx1, idx2) = min(tuple([dist, idx1, idx2]), closest_pair_strip(cluster_list, mid, dist),
                                 key=lambda tup: tup[0])
    return tuple([dist, idx1, idx2])


def hierarchical_clustering(cluster_list, num_clusters):
    """creating clusters using hierarchical structure"""
    clusters = [cluster for cluster in cluster_list]
    while len(clusters) > num_clusters:
        clstr_lst = list(clusters)
        clstr_lst.sort(key=lambda cls: cls.horiz_center())
        closest = fast_closest_pair(clstr_lst)
        clusters.append(clstr_lst[closest[1]].merge_clusters(clstr_lst[closest[2]]))
        clusters.remove(clstr_lst[closest[1]])
        clusters.remove(clstr_lst[closest[2]])
    return clusters


def kmeans_clustering(cluster_list, num_clusters, num_iterations):
    """implementation of kmeans clustering"""
    clusters = [cluster.copy() for cluster in cluster_list]
    clusters.sort(key=lambda cls: cls.total_population(), reverse=True)
    centers = [center for center in clusters[:num_clusters]]
    for dummy_iteration in range(num_iterations):
        temps = [[] for idx in range(num_clusters)]
        for idx in range(len(cluster_list)):
            dist = float("Inf")
            for center in range(num_clusters):
                if cluster_list[idx].distance(centers[center]) < dist:
                    dist = cluster_list[idx].distance(centers[center])
                    cluster_sel = center
            if not temps[cluster_sel]:
                temps[cluster_sel] = cluster_list[idx].copy()
            else:
                temps[cluster_sel].merge_clusters(cluster_list[idx])
        for center in range(num_clusters):
            centers[center] = temps[center].copy()
    return temps


# cluster_lst = [Cluster(set(['01073']), 704.191210749, 411.014665198, 662047, 7.3e-05), Cluster(set(['06059']),
# 113.997715586, 368.503452566, 2846289, 9.8e-05), Cluster(set(['06037']), 105.369854549, 359.050126004, 9519338,
# 0.00011), Cluster(set(['06029']), 103.787886113, 326.006585349, 661645, 9.7e-05), Cluster(set(['06071']),
# 148.402461892, 350.061039619, 1709434, 7.7e-05), Cluster(set(['06075']), 52.7404001225, 254.517429395, 776733,
# 8.4e-05), Cluster(set(['08031']), 371.038986573, 266.847932979, 554636, 7.9e-05), Cluster(set(['24510']),
# 872.946822486, 249.834427518, 651154, 7.4e-05), Cluster(set(['34013']), 906.236730753, 206.977429459, 793633,
# 7.1e-05), Cluster(set(['34039']), 905.587082153, 210.045085725, 522541, 7.3e-05), Cluster(set(['34017']),
# 909.08042421, 207.462937763, 608975, 9.1e-05), Cluster(set(['36061']), 911.072622034, 205.783086757, 1537195,
# 0.00015), Cluster(set(['36005']), 912.315497328, 203.674106811, 1332650, 0.00011), Cluster(set(['36047']),
# 911.595580089, 208.928374072, 2465326, 9.8e-05), Cluster(set(['36059']), 917.384980291, 205.43647538, 1334544,
# 7.6e-05), Cluster(set(['36081']), 913.462051588, 207.615750359, 2229379, 8.9e-05), Cluster(set(['41051']),
# 103.293707198, 79.5194104381, 660486, 9.3e-05), Cluster(set(['41067']), 92.2254623376, 76.2593957841, 445342,
# 7.3e-05), Cluster(set(['51013']), 865.681962839, 261.222875114, 189453, 7.7e-05), Cluster(set(['51840']),
# 845.843602685, 258.214178983, 23585, 7.1e-05), Cluster(set(['51760']), 865.424050159, 293.735963553, 197790,
# 8.6e-05), Cluster(set(['55079']), 664.855000617, 192.484141264, 940164, 7.4e-05), Cluster(set(['54009']),
# 799.221537984, 240.153315109, 25447, 7.7e-05), Cluster(set(['11001']), 867.470401202, 260.460974222, 572059,
# 7.7e-05)] print hierarchical_clustering(cluster_lst, 23)
