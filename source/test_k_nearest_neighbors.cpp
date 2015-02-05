#include <macros.hpp>
#include <k_nearest_neighbors.hpp>
#include <numerics.hpp>
#include <test_k_nearest_neighbors.hpp>

std::string
TestKNearestNeighbors::name() const
{
    return std::string("K NEAREST NEIGHBORS TESTER");
}

unsigned int
TestKNearestNeighbors::numberOfSubtests() const
{
    return 1;
}

void
TestKNearestNeighbors::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < 1, "invalid index " << i);

    std::vector<std::vector<double> > points;

    for(int i = -100; i < 100; ++i)
    {
        for(int j = -100; j < 100; ++j)
        {
            int n = points.size();
            points.resize(n + 1);

            points[n].resize(2);
            points[n][0] = double(i);
            points[n][1] = double(j);
        }
    }

    Math::KNearestNeighbors knn(2, points, 3);

    std::vector<double> q(2);
    q[0] = 50.1;
    q[1] = 20.2;

    std::vector<unsigned long> indices;
    std::vector<double> distances;

    knn.search(q, &indices, &distances);

    res = 1;
    expected = 1;
    subTestName = "2d_grid";

    if(points[indices[0]][0] != 50 || points[indices[0]][1] != 20)
    {
        output_screen("FAIL: First neighbor should be (50, 20) but it is (" << points[indices[0]][0] << ", " << points[indices[0]][1] << ")." << std::endl);
        res = 0;
    }

    if(points[indices[1]][0] != 50 || points[indices[1]][1] != 21)
    {
        output_screen("FAIL: Second neighbor should be (50, 21) but it is (" << points[indices[1]][0] << ", " << points[indices[1]][1] << ")." << std::endl);
        res = 0;
    }

    if(points[indices[2]][0] != 51 || points[indices[2]][1] != 20)
    {
        output_screen("FAIL: Third neighbor should be (51, 20) but it is (" << points[indices[2]][0] << ", " << points[indices[2]][1] << ")." << std::endl);
        res = 0;
    }

    if(!Math::areEqual(distances[0], 0.05, 1e-7))
    {
        output_screen("FAIL: Distance squared to the first neighbor should be " << 0.05 << " but it is " << distances[0] << "." << std::endl);
        res = 0;
    }
}
