#ifndef COSMO_PP_KD_TREE_HPP
#define COSMO_PP_KD_TREE_HPP

#include <vector>
#include <queue>
#include <utility>

/// A k-d tree class.
class KDTree
{
public:
    /// Constructor.
    /// \param dim The dimensionality of the space.
    /// \param points The points to build the tree on.
    KDTree(int dim, const std::vector<std::vector<double> >& points);

    /// Destructor.
    ~KDTree();

    /// Reset the set of points and rebuild the tree.
    void reset(const std::vector<std::vector<double> >& points);

    /// Rebalance the tree. The set of points doesn't change.
    void reBalance();

    /// Insert a new point into the tree. The tree is not re-balanced, the element is simply inserted into the appropriate location of the tree.
    /// \param point The point to insert.
    void insert(const std::vector<double>& point);

    /// Get the depth of the tree.
    /// \return The depth.
    int depth() const { return depth_; }

    /// Get the number of elements in the tree.
    /// \return The number of elements.
    unsigned long nElements() const { return elements_.size(); }

    /// Find k nearest neighbors of a given point in the tree.
    /// \param point The point to search around.
    /// \param k The number of nearest neighbors to return.
    /// \param neighbors A pointer to a vector where the nearest neighbors will be returned.
    /// \param distanceSquares A pointer to a vector where the squared distances to the nearest neighbors will be returned. Can be set to NULL (default option) in which case this will be ignored.
    void findNearestNeighbors(const std::vector<double> &point, int k, std::vector<std::vector<double> > *neighbors, std::vector<double> *distanceSquares = NULL) const;

    /// Find k nearest neighbors of a given point in the tree.
    /// \param point The point to search around.
    /// \param k The number of nearest neighbors to return.
    /// \indices A pointer to a vector where the INDICES of the nearest neighbors will be returned, as they were in the original vector in the constructor. Note that the index of an element added by the insert function will simply follow the highest index already in the tree.
    /// \param distanceSquares A pointer to a vector where the squared distances to the nearest neighbors will be returned. Can be set to NULL (default option) in which case this will be ignored.
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

