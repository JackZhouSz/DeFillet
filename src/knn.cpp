//
// Created by xiaowuga on 12/16/2023.
//

#include "knn.h"

namespace KNN {

// 定义 KD 树适配器
    class PointCloud {
    public:
        PointCloud(const std::vector<Point>* points) : pts(points) {}
        const std::vector<Point>* pts;

        // 提取点云数据
        inline size_t kdtree_get_point_count() const { return pts->size(); }

        // 访问第 idx 个点的坐标
        inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
            return pts->at(idx).p[dim];
        }

        // 回调函数：计算点之间的距离的平方
        template<class BBOX>
        bool kdtree_get_bbox(BBOX &) const { return false; }
    };

    struct KdTree : public nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, PointCloud>, PointCloud, 3, int> {
        KdTree(PointCloud* cloud) :
        KDTreeSingleIndexAdaptor< nanoflann::L2_Simple_Adaptor<float, PointCloud>, PointCloud, 3, int >(3,
                                                                                                         *cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10)){}
        ~KdTree() {
//            delete cloud_;
        }

        PointCloud* cloud_;
    };

    #define get_tree(x) (reinterpret_cast<const KdTree *>(x))

    KdSearch::KdSearch(std::vector<Point>& points) {
        points_ =  const_cast< std::vector<Point>* >(&points);

        PointCloud* cloud = new PointCloud(points_);
        KdTree* tree = new KdTree(cloud);
        tree->buildIndex();
        tree_ = tree;
    }

    KdSearch::~KdSearch() {
        delete get_tree(tree_);
    }

    int KdSearch::radius_search(const Point& p, double radius, std::vector<size_t> &neighbors,
                                 std::vector<float> &squared_distances) const {
        std::vector<std::pair<int , float> > matches;
        nanoflann::SearchParams params;
        params.sorted = true;
        const std::size_t num = get_tree(tree_)->radiusSearch(p.p, radius, matches, params);

        neighbors.resize(num);
        squared_distances.resize(num);
        for (std::size_t i = 0; i < num; ++i) {
            const std::pair<std::size_t, float>& e = matches[i];
            neighbors[i] = e.first;
            squared_distances[i] = e.second;
        }

        return (int)neighbors.size();
    }

    void KdSearch::kth_search(const Point &p, int k, std::vector<size_t> &neighbors,
                              std::vector<float> &squared_distances) const {
        std::vector<int> indices(k);
        std::vector<float> distances(k);

        // float distances[k];
        // result_set.init(&indices[0], &sqr_distances[0]);
        get_tree(tree_)->knnSearch(&p.p[0], k, &indices[0], &distances[0]);
        neighbors.resize(k);
        squared_distances.resize(k);
        for(int i = 0; i < k; i++) {
            neighbors[i] = indices[i];
            squared_distances[i] = distances[i];
        }
        // get_tree(tree_)->findNeighbors(result_set, p.p,  nanoflann::SearchParams(k));
    }
}