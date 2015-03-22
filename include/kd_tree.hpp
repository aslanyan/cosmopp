#ifndef COSMO_PP_KD_TREE_HPP
#define COSMO_PP_KD_TREE_HPP

#include <vector>
#include <queue>
#include <utility>

class KDTree
{
public:
    KDTree(int dim, const std::vector<std::vector<double> >& elements);

    ~KDTree();

    void reset(const std::vector<std::vector<double> >& points);
    void reBalance();

    void insert(const std::vector<double>& elem);

    int depth() const { return depth_; }
    unsigned long nElements() const { return elements_.size(); }

    void findNearestNeighbors(const std::vector<double> &point, int k, std::vector<std::vector<double> > *neighbors, std::vector<double> *distanceSquares = NULL) const;
    void findNearestNeighbors(const std::vector<double> &point, int k, std::vector<unsigned long> *indices, std::vector<double> *distanceSquares = NULL) const;

private:
    struct Node
    {
        unsigned long index;
        int depth;
        Node *parent;
        Node *left;
        Node *right;
    };

    struct ComparePair
    {
        inline bool operator() (const std::pair<double, unsigned long>& a, const std::pair<double, unsigned long>& b) const
        {
            return a.first < b.first;
        }
    };

private:
    Node* construct(std::vector<unsigned long> &elementsIndices, unsigned long begin, unsigned long end, int depth, Node *parent);
    void destruct(Node *root);

    void search(const Node *current, std::priority_queue<std::pair<double, unsigned long>, std::vector<std::pair<double, unsigned long> >, ComparePair>& bpq, int k, const std::vector<double>& point) const;

private:
    int dim_;
    std::vector<std::vector<double> > elements_;

    int depth_;

    Node *root_;
};

#endif

