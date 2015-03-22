#include <algorithm>

#include <macros.hpp>
#include <kd_tree.hpp>

KDTree::KDTree(int dim, const std::vector<std::vector<double> >& elements) : dim_(dim), elements_(elements)
{
    check(dim_ > 0, "invalid dimension " << dim_ << ", must be positive");

    std::vector<unsigned long> elemsIndices(elements_.size());
    for(unsigned long i = 0; i < elemsIndices.size(); ++i)
        elemsIndices[i] = i;

    depth_ = 0;
    root_ = construct(elemsIndices, 0, elemsIndices.size(), 0, NULL);
}

KDTree::~KDTree()
{
    destruct(root_);
}

void
KDTree::reset(const std::vector<std::vector<double> >& elements)
{
    elements_ = elements;
    reBalance(); // assumes that rebalance destroys and reconstructs the tree
}

void
KDTree::reBalance()
{
    destruct(root_);

    std::vector<unsigned long> elemsIndices(elements_.size());
    for(unsigned long i = 0; i < elemsIndices.size(); ++i)
        elemsIndices[i] = i;

    depth_ = 0;
    root_ = construct(elemsIndices, 0, elemsIndices.size(), 0, NULL);
}

void
KDTree::insert(const std::vector<double>& elem)
{
    check(elem.size() == dim_, "");

    elements_.push_back(elem);

    Node *node = new Node;
    node->index = elements_.size() - 1;
    node->left = NULL;
    node->right = NULL;

    if(!root_)
    {
        check(elements_.size() == 1, "");
        node->depth = 0;
        node->parent = NULL;
        root_ = node;

        check(depth_ == 0, "");
        depth_ = 1;
        return;
    }

    Node *current = root_;

    while(true)
    {
        const int compareIndex = current->depth % dim_;
        const bool goLeft = elem[compareIndex] < elements_[current->index][compareIndex];

        if(goLeft)
        {
            if(current->left)
            {
                current = current->left;
                continue;
            }
            else
            {
                current->left = node;
                break;
            }
        }
        
        check(!goLeft, "");
        if(current->right)
        {
            current = current->right;
            continue;
        }

        check(!current->right, "");
        current->right = node;
        break;
    }

    node->parent = current;
    node->depth = current->depth + 1;
}

void
KDTree::findNearestNeighbors(const std::vector<double> &point, int k, std::vector<std::vector<double> > *neighbors, std::vector<double> *distanceSquares) const
{
    std::vector<unsigned long> indices;
    findNearestNeighbors(point, k, &indices, distanceSquares);

    check(indices.size() == k, "");

    neighbors->resize(k);
    for(int i = 0; i < k; ++i)
    {
        check(indices[i] < elements_.size(), "");
        neighbors->at(i) = elements_[indices[i]];
    }
}

void
KDTree::findNearestNeighbors(const std::vector<double>& point, int k, std::vector<unsigned long> *indices, std::vector<double> *distanceSquares) const
{
    check(point.size() == dim_, "");

    check(k >= 0, "invalid k");
    check(indices, "");

    indices->resize(k);

    if(distanceSquares)
        distanceSquares->resize(k);

    if(!k)
        return;

    check(nElements() >= k, k << " nearest neighbors requested but there are only " << nElements() << " elements in the kd tree");

    ComparePair cp;
    std::vector<std::pair<double, unsigned long> > container;
    container.reserve(k + 1);
    std::priority_queue<std::pair<double, unsigned long>, std::vector<std::pair<double, unsigned long> >, ComparePair> bpq(cp, container);
    search(root_, bpq, k, point);

    check(bpq.size() == k, "");

    for(int i = k - 1; i >= 0; --i)
    {
        check(!bpq.empty(), "");
        indices->at(i) = bpq.top().second;
        if(distanceSquares)
            distanceSquares->at(i) = bpq.top().first;
        bpq.pop();
    }

    check(bpq.empty(), "");
}

void
KDTree::search(const Node *current, std::priority_queue<std::pair<double, unsigned long>, std::vector<std::pair<double, unsigned long> >, ComparePair>& bpq, int k, const std::vector<double>& point) const
{
    check(k > 0, "");
    check(point.size() == dim_, "");

    if(!current)
        return;

    check(current->index < elements_.size(), "");
    const std::vector<double>& v = elements_[current->index];
    check(v.size() == dim_, "");

    // calculate distance
    double distance = 0;
    for(int i = 0; i < dim_; ++i)
    {
        const double delta = point[i] - v[i];
        distance += delta * delta; // euclidean distance squared
    }

    check(bpq.size() <= k, "");

    bpq.emplace(distance, current->index);
    if(bpq.size() > k)
        bpq.pop();

    check(bpq.size() <= k, "");

    const int index = current->depth % dim_;
    const bool goLeft = point[index] < v[index];

    if(goLeft)
        search(current->left, bpq, k, point);
    else
        search(current->right, bpq, k, point);

    bool goOtherSide = false;
    if(bpq.size() < k)
        goOtherSide = true;
    else
    {
        const double delta = point[index] - v[index];
        const double deltaSq = delta * delta;

        check(!bpq.empty(), "");
        const double maxDist = bpq.top().first;

        if(deltaSq <= maxDist)
            goOtherSide = true;
    }

    if(!goOtherSide)
        return;

    // search the other branch
    if(goLeft)
        search(current->right, bpq, k, point);
    else
        search(current->left, bpq, k, point);
}

namespace
{

struct CompareKDTreeNode
{
    int dim;
    int compareIndex;
    const std::vector<std::vector<double> > *elements;

    bool operator() (unsigned long a, unsigned long b) const
    {
        check(dim > 0, "");
        check(compareIndex < dim, "");

        check(a < elements->size(), "");
        check(b < elements->size(), "");

        const std::vector<double> &aVec = elements->at(a);
        const std::vector<double> &bVec = elements->at(b);

        check(aVec.size() == dim, "");
        check(bVec.size() == dim, "");

        return aVec[compareIndex] < bVec[compareIndex];
    }
};

}

KDTree::Node*
KDTree::construct(std::vector<unsigned long> &elementsIndices, unsigned long begin, unsigned long end, int depth, Node *parent)
{
    check(end >= begin, "");
    check(end <= elementsIndices.size(), "");

    if(end == begin)
        return NULL;

    CompareKDTreeNode comp;
    comp.dim = dim_;
    comp.compareIndex = depth % dim_;
    comp.elements = &elements_;

    std::vector<unsigned long>::iterator beginIt = elementsIndices.begin() + begin, endIt = elementsIndices.begin() + end;
    std::sort(beginIt, endIt, comp);

    unsigned long median = (begin + end) / 2;
    check(median >= begin && median < end, "");

    Node* node = new Node;
    node->index = elementsIndices[median];
    node->depth = depth;
    node->parent = parent;
    node->left = construct(elementsIndices, begin, median, depth + 1, node);
    node->right = construct(elementsIndices, median + 1, end, depth + 1, node);

    depth_ = std::max(depth + 1, depth_);

    return node;
}

void
KDTree::destruct(Node* root)
{
    if(!root)
        return;

    destruct(root->left);
    destruct(root->right);
    delete root;
}

