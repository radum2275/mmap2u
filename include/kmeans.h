/*
 * potential.h
 *
 *  Created on: Apr 12, 2023
 *      Author: radu
 *
 * Copyright (c) 2015, International Business Machines Corporation
 * and University of California Irvine. All rights reserved.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


/// \file potential.h
/// \brief A table based potential for credal networks
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_LOOPY_KMEANS_H_
#define IBM_LOOPY_KMEANS_H_

#include "base.h"

namespace merlin {

class point {
protected:
    size_t m_id; 
    size_t m_cluster;
    size_t m_dimensions;
    std::vector<double> m_values;

public:
    point(size_t id, const std::vector<double>& values) : m_id(id), m_values(values) {
        m_dimensions = values.size();
        m_cluster = 0; // Initially not assigned to any cluster
    }
    point(const point& p) {
        m_id = p.m_id;
        m_cluster = p.m_cluster;
        m_dimensions = p.m_dimensions;
        m_values = p.m_values;
    }
    point& operator=(const point& p) {
        m_id = p.m_id;
        m_cluster = p.m_cluster;
        m_dimensions = p.m_dimensions;
        m_values = p.m_values;
        return *this;
    }

    inline size_t get_dimensions() { return m_dimensions; }
    inline size_t get_cluster() { return m_cluster; }
    inline size_t get_id() { return m_id; }
    inline void set_cluster(size_t val) { m_cluster = val; }
    inline double get_val(size_t pos) { return m_values[pos]; }
};
    
class cluster {
protected:
    int m_id;
    std::vector<double> m_centroid;
    std::vector<point> m_points;

public:
    cluster(int id, point& centroid) {
        m_id = id;
        for (size_t i = 0; i < centroid.get_dimensions(); ++i) {
            m_centroid.push_back(centroid.get_val(i));
        }
        this->add_point(centroid);
    }

    void add_point(point& p) {
        p.set_cluster(this->m_id);
        m_points.push_back(p);
    }

    bool remove_point(size_t pid) {
        size_t num_points = m_points.size();

        for (size_t i = 0; i < num_points; i++) {
            if (m_points[i].get_id() == pid) {
                m_points.erase(m_points.begin() + i);
                return true;
            }
        }

        return false;
    }

    void remove_all() { m_points.clear(); }

    int get_id() { return m_id; }

    point& get_point(size_t pos) { return m_points[pos]; }

    size_t size() { return m_points.size(); }

    double get_centroid_by_pos(size_t pos) { return m_centroid[pos]; }

    void set_centroid_by_pos(size_t pos, double val) { this->m_centroid[pos] = val; }
};

class kmeans {
protected:
    size_t m_k;
    size_t m_iterations;
    size_t m_dimensions; 
    size_t m_total_points;
    std::string m_distance;
    std::vector<cluster> m_clusters;
    std::vector<std::vector<double> > m_bounds;

    void clear_clusters() {
        for (size_t i = 0; i < m_k; i++) {
            m_clusters[i].remove_all();
        }
    }

    size_t get_nearest_cluster(point& p) {
        double sum = 0.0, min_dist;
        size_t nearest_cluster;
        if (m_dimensions == 1) {
            min_dist = std::abs(m_clusters[0].get_centroid_by_pos(0) - p.get_val(0));
        } else  {
            for (size_t i = 0; i < m_dimensions; i++) {
                if (m_distance == "euclidian") {
                    sum += std::pow(m_clusters[0].get_centroid_by_pos(i) - p.get_val(i), 2.0);
                } else if (m_distance == "manhattan") {
                    sum += std::abs(m_clusters[0].get_centroid_by_pos(i) - p.get_val(i));
                }
            }
            
            if (m_distance == "euclidian") {
                min_dist = sqrt(sum);
            } else if (m_distance == "manhattan") {
                min_dist = sum;
            }
        }
        
        nearest_cluster = m_clusters[0].get_id();

        for (size_t i = 1; i < m_k; i++) {
            double dist;
            sum = 0.0;
            
            if (m_dimensions == 1) {
                dist = std::abs(m_clusters[i].get_centroid_by_pos(0) - p.get_val(0));
            } else {
                for (size_t j = 0; j < m_dimensions; j++) {
                    if (m_distance == "euclidian") {
                        sum += std::pow(m_clusters[i].get_centroid_by_pos(j) - p.get_val(j), 2.0);
                    } else if (m_distance == "manhattan") {
                        sum += std::abs(m_clusters[i].get_centroid_by_pos(j) - p.get_val(j));
                    }
                }

                if (m_distance == "euclidian") {
                    dist = sqrt(sum);
                } else if (m_distance == "manhattan") {
                    dist = sum;
                }
            }

            if (dist < min_dist) {
                min_dist = dist;
                nearest_cluster = m_clusters[i].get_id();
            }
        }

        return nearest_cluster;
    }

public:
    
    kmeans(size_t k, size_t iterations, std::string distance) {
        m_k = k;
        m_iterations = iterations;
        m_distance = distance;
        srand(12345678);
    }

    std::vector<std::vector<double> >& get_bounds() {
        return m_bounds;
    }

    void run(std::vector<point> &all_points) {
        m_total_points = all_points.size();
        m_dimensions = all_points[0].get_dimensions();

        // Initializing Clusters
        std::vector<size_t> used_point_ids;

        for (size_t i = 1; i <= m_k; i++) {
            while (true) {
                int index = rand() % m_total_points;

                if (find(used_point_ids.begin(), used_point_ids.end(), index) == used_point_ids.end()) {
                    used_point_ids.push_back(index);
                    all_points[index].set_cluster(i);
                    cluster cl(i, all_points[index]);
                    m_clusters.push_back(cl);
                    break;
                }
            }
        }

        size_t iter = 1;
        while (true) {
            bool done = true;

            // Add all points to their nearest cluster
            for (size_t i = 0; i < m_total_points; i++) {
                size_t current_cluster = all_points[i].get_cluster();
                size_t nearest_cluster = get_nearest_cluster(all_points[i]);

                if (current_cluster != nearest_cluster) {
                    all_points[i].set_cluster(nearest_cluster);
                    done = false;
                }
            }

            // clear all existing clusters
            clear_clusters();

            // reassign points to their new clusters
            for (size_t i = 0; i < m_total_points; i++) {
                // cluster index is ID-1
                m_clusters[all_points[i].get_cluster() - 1].add_point(all_points[i]);
            }

            // Recalculating the center of each cluster
            for (size_t i = 0; i < m_k; i++) {
                size_t cluster_size = m_clusters[i].size();

                for (size_t j = 0; j < m_dimensions; j++) {
                    double sum = 0.0;
                    if (cluster_size > 0) {
                        for (size_t p = 0; p < cluster_size; p++) {
                            sum += m_clusters[i].get_point(p).get_val(j);
                        }
                        m_clusters[i].set_centroid_by_pos(j, sum / cluster_size);
                    }
                }
            }

            if (done || iter >= m_iterations) {
                break;
            }

            iter++;
        }

        // For each cluster, get the Pareto Leaset Upper Bound (max)
        for (size_t i = 0; i < m_k; i++) {
            std::vector<double> plub(m_dimensions, 0.0);
            for (size_t p = 0; p < m_clusters[i].size(); ++p) {
                for (size_t j = 0; j < m_dimensions; ++j) {
                    plub[j] = std::max(plub[j], m_clusters[i].get_point(p).get_val(j));
                }
            }
            m_bounds.push_back(plub);
        }
    }
};



} // end namespace


#endif
