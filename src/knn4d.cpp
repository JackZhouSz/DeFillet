//
// Created by xiaowuga on 2025/8/21.
//

#include <knn4d.h>

namespace KNN {

// 定义 KD 树适配器
    class PointCloud4D {
    public:
        PointCloud4D(const std::vector<Point4D>* points) : pts(points) {}
        const std::vector<Point4D>* pts;

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

    struct KdTree4D : public nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud4D>, PointCloud4D, 4, int> {
        KdTree4D(PointCloud4D* cloud) :
        KDTreeSingleIndexAdaptor< nanoflann::L2_Simple_Adaptor<double, PointCloud4D>, PointCloud4D, 4, int >(4,
                                                                                                         *cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10)){}
        ~KdTree4D() {
//            delete cloud_;
        }

        PointCloud4D* cloud_;
    };

    #define get_tree4D(x) (reinterpret_cast<const KdTree4D *>(x))

    KdSearch4D::KdSearch4D(std::vector<Point4D>& points) {
        points_ =  const_cast< std::vector<Point4D>* >(&points);

        PointCloud4D* cloud = new PointCloud4D(points_);
        KdTree4D* tree = new KdTree4D(cloud);
        tree->buildIndex();
        tree_ = tree;
    }

    KdSearch4D::~KdSearch4D() {
        delete get_tree4D(tree_);
    }

    // int KdSearch4D::radius_search(const Point4D& p, double radius, std::vector<size_t> &neighbors,
    //                              std::vector<double> &squared_distances) const {
    //     std::vector<std::pair<int , double> > matches;
    //     nanoflann::SearchParams params;
    //     params.sorted = true;
    //     const std::size_t num = get_tree(tree_)->radiusSearch(p.p, radius, matches, params);
    //
    //     neighbors.resize(num);
    //     squared_distances.resize(num);
    //     for (std::size_t i = 0; i < num; ++i) {
    //         const std::pair<std::size_t, float>& e = matches[i];
    //         neighbors[i] = e.first;
    //         squared_distances[i] = e.second;
    //     }
    //
    //     return (int)neighbors.size();
    // }

    void KdSearch4D::kth_search(const Point4D &p, int k, std::vector<size_t> &neighbors,
                              std::vector<double> &squared_distances) const {
        std::vector<int> indices(k);
        std::vector<double> distances(k);

        // float distances[k];
        // result_set.init(&indices[0], &sqr_distances[0]);
        get_tree4D(tree_)->knnSearch(&p.p[0], k, &indices[0], &distances[0]);
        neighbors.resize(k);
        squared_distances.resize(k);
        for(int i = 0; i < k; i++) {
            neighbors[i] = indices[i];
            squared_distances[i] = distances[i];
        }
        // get_tree(tree_)->findNeighbors(result_set, p.p,  nanoflann::SearchParams(k));
    }
}