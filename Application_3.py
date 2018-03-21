"""Analysis of cancer risk"""
import math
import random
import time
import csv
import urllib2
from matplotlib import pyplot as plt


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


def scale_random():
    num = random.random()
    if num < 0.5:
        return -2 * num
    else:
        return (num - 0.5) * 2


def gen_random_clusters(num_clusters):
    clusters = []
    for _ in range(num_clusters):
        clusters.append(Cluster('', scale_random(), scale_random(), 0, 0))

    return clusters


def compare_running_times():
    time1 = time.time()
    slow_times = []
    fast_times = []

    for size in range(2, 201):
        cluster = gen_random_clusters(size)

        start = time.time()
        slow_closest_pair(cluster)
        stop = time.time()
        slow_times.append(stop - start)

        start = time.time()
        fast_closest_pair(cluster)
        stop = time.time()
        fast_times.append(stop - start)

        print "Done with node of size ", size

    plt.plot([n for n in range(2, 201)], slow_times, 'r', label="slow_closest_pair")
    plt.plot([n for n in range(2, 201)], fast_times, 'b', label="fast_closest_pair")
    plt.legend(loc="upper left")
    plt.ylabel("Running Time (sec)")
    plt.xlabel("Number of Nodes in Input")
    plt.title("Comparison of Running Times of Closest Pair Algorithms")
    time2 = time.time()
    print "Total took ", time2 - time1, " seconds"
    plt.show()


# URLS for various important datasets
DIRECTORY = "http://commondatastorage.googleapis.com/codeskulptor-assets/"
MAP_URL = DIRECTORY + "data_clustering/USA_Counties.png"

# Define colors for clusters.  Display a max of 16 clusters.
COLORS = ['Aqua', 'Yellow', 'Blue', 'Fuchsia', 'Black', 'Green', 'Lime', 'Maroon', 'Navy', 'Olive', 'Orange', 'Purple',
          'Red', 'Brown', 'Teal']


# Helper functions

def circle_area(pop):
    """
    Compute area of circle proportional to population
    """
    return math.pi * pop / (200.0 ** 2)


def plot_clusters(data_table, cluster_list, draw_centers=False):
    """
    Create a plot of clusters of counties
    """

    fips_to_line = {}
    for line_idx in range(len(data_table)):
        fips_to_line[data_table[line_idx][0]] = line_idx

    # Load map image
    map_file = urllib2.urlopen(MAP_URL)
    map_img = plt.imread(map_file)

    # Scale plot to get size similar to CodeSkulptor version
    ypixels, xpixels, bands = map_img.shape
    DPI = 60.0  # adjust this constant to resize your plot
    xinch = xpixels / DPI
    yinch = ypixels / DPI
    plt.figure(figsize=(xinch, yinch))
    implot = plt.imshow(map_img)

    # draw the counties colored by cluster on the map
    if not draw_centers:
        for cluster_idx in range(len(cluster_list)):
            cluster = cluster_list[cluster_idx]
            cluster_color = COLORS[cluster_idx % len(COLORS)]
            for fips_code in cluster.fips_codes():
                line = data_table[fips_to_line[fips_code]]
                plt.scatter(x=[line[1]], y=[line[2]], s=circle_area(line[3]), lw=1,
                            facecolors=cluster_color, edgecolors=cluster_color)

    # add cluster centers and lines from center to counties
    else:
        for cluster_idx in range(len(cluster_list)):
            cluster = cluster_list[cluster_idx]
            cluster_color = COLORS[cluster_idx % len(COLORS)]
            for fips_code in cluster.fips_codes():
                line = data_table[fips_to_line[fips_code]]
                plt.scatter(x=[line[1]], y=[line[2]], s=circle_area(line[3]), lw=1,
                            facecolors=cluster_color, edgecolors=cluster_color, zorder=1)
        for cluster_idx in range(len(cluster_list)):
            cluster = cluster_list[cluster_idx]
            cluster_color = COLORS[cluster_idx % len(COLORS)]
            cluster_center = (cluster.horiz_center(), cluster.vert_center())
            for fips_code in cluster.fips_codes():
                line = data_table[fips_to_line[fips_code]]
                plt.plot([cluster_center[0], line[1]], [cluster_center[1], line[2]], cluster_color, lw=1, zorder=2)
        for cluster_idx in range(len(cluster_list)):
            cluster = cluster_list[cluster_idx]
            cluster_color = COLORS[cluster_idx % len(COLORS)]
            cluster_center = (cluster.horiz_center(), cluster.vert_center())
            cluster_pop = cluster.total_population()
            plt.scatter(x=[cluster_center[0]], y=[cluster_center[1]], s=circle_area(cluster_pop), lw=2,
                        facecolors="none", edgecolors="black", zorder=3)

    plt.show()


def create_data_from_csv(file_name):
    result = []
    with open(file_name, 'rb') as file:
        reader = csv.reader(file)
        for row in reader:
            params = [item for item in row]
            result.append([params[0], float(params[1]), float(params[2]), int(params[3]), float(params[4])])
            # params = [item for item in row]
            # result.append(Cluster(set(params[0]), float(params[1]), float(params[2]), int(params[3]), float(params[4])))

    return result


def create_clusters_from_data(data):
    result = []
    for param in data:
        result.append(Cluster(set(param[0]), float(param[1]), float(param[2]), int(param[3]), float(param[4])))

    return result


DIRECTORY = "http://commondatastorage.googleapis.com/codeskulptor-assets/"
DATA_3108_URL = DIRECTORY + "data_clustering/unifiedCancerData_3108.csv"
DATA_896_URL = DIRECTORY + "data_clustering/unifiedCancerData_896.csv"
DATA_290_URL = DIRECTORY + "data_clustering/unifiedCancerData_290.csv"
DATA_111_URL = DIRECTORY + "data_clustering/unifiedCancerData_111.csv"


def load_data_table(data_url):
    """
    Import a table of county-based cancer risk data
    from a csv format file
    """
    data_file = urllib2.urlopen(data_url)
    data = data_file.read()
    data_lines = data.split('\n')
    print "Loaded", len(data_lines), "data points"
    data_tokens = [line.split(',') for line in data_lines]
    return [[tokens[0], float(tokens[1]), float(tokens[2]), int(tokens[3]), float(tokens[4])]
            for tokens in data_tokens]


# data = load_data_table(DATA_111_URL)
# print "Done creating data"
# singleton_list = [Cluster({line[0]}, line[1], line[2], line[3], line[4]) for line in data]
# hclusters = hierarchical_clustering(singleton_list, 9)
# kclusters = kmeans_clustering(singleton_list, 9, 5)
# print "Done clustering"
# plot_clusters(data, clusters, True)


def compute_distortion(cluster_list, data_table):
    distortion = 0
    for cluster in cluster_list:
        distortion += cluster.cluster_error(data_table)

    return distortion

# print compute_distortion(hclusters, data)
# print compute_distortion(kclusters, data)


def compare_distortions(data_table):
    data = load_data_table(data_table)
    singletons = [Cluster({line[0]}, line[1], line[2], line[3], line[4]) for line in data]
    kdistortions = []
    hdistortions = []
    for num in range(6, 21):
        singletons = [Cluster({line[0]}, line[1], line[2], line[3], line[4]) for line in data]
        kclustering = kmeans_clustering(singletons, num, 5)
        kdistortions.append(compute_distortion(kclustering, data))

        singletons = [Cluster({line[0]}, line[1], line[2], line[3], line[4]) for line in data]
        hclustering = hierarchical_clustering(singletons, num)
        hdistortions.append(compute_distortion(hclustering, data))

        print "Done with iteration ", num

    plt.plot([num for num in range(6, 21)], kdistortions, 'r', label="k-means clustering")
    plt.plot([num for num in range(6, 21)], hdistortions, 'b', label="hierarchical clustering")
    plt.legend(loc="upper right")
    plt.xlabel("Number of Output Nodes")
    plt.ylabel("Distortion")
    plt.title("Distortion of Clustering Algorithms - 896 Counties")
    plt.show()

compare_distortions(DATA_896_URL)

