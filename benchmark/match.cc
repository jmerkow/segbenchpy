
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include "csa.hh"
#include "kofn.hh"
#include "Point.hh"
#include "Matrix.hh"
#include "Array.hh"
#include "match.hh"
#include "Timer.hh"

struct Edge {
    int i,j;	// node ids, 0-based
    double w;	// distance between pixels
};

// CSA code needs integer weights.  Use this multiplier to convert
// floating-point weights to integers.
static const int multiplier = 1000;

// The degree of outlier connections.
// static const int degree = 6;



double matchEdgeMaps2D(const double* bmap1, const int height, const int width, 
                     const double* bmap2, const int height2, const int width2,
                     double* m1, int m1h, int m1w, 
                     double* m2, int m2h, int m2w,
                     double maxDist, double outlierCost, int degree)
{
    // Check global constants.
    assert (degree > 0);
    assert (multiplier > 0);

    // Check arguments.
    assert (height == height2);
    assert (width == width2);
    assert (maxDist >= 0);
    assert (outlierCost > maxDist);

    // Initialize match[12] arrays to (-1,-1).
    std::vector< std::vector< Pixel > > match1 (height , std::vector< Pixel > ( width, Pixel(-1,-1) ));
    std::vector< std::vector< Pixel > > match2 (height , std::vector< Pixel > ( width, Pixel(-1,-1) ));
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            match1[y][x] = Pixel(-1,-1);
            match2[y][x] = Pixel(-1,-1);
        }
    }

    // Radius of search window.
    const int r = (int) ceil (maxDist); 

    // Figure out which nodes are matchable, i.e. within maxDist
    // of another node.
    std::vector< std::vector< int > > matchable1 (height, std::vector< int > ( width, 0 ));
    std::vector< std::vector< int > > matchable2 (height, std::vector< int > ( width, 0 ));


    for (int y1 = 0; y1 < height; ++y1) {
        for (int x1 = 0; x1 < width; ++x1) {
            if (!bmap1[y1*width+x1]) { continue; }
            for (int v = -r; v <= r; ++v) {
                for (int u = -r; u <= r; ++u) {
                    const double d2 = u*u + v*v;
                    if (d2 > maxDist*maxDist) { continue; }
                    const int x2 = x1 + u;
                    const int y2 = y1 + v;
                    if (x2 < 0 || x2 >= width) { continue; }
                    if (y2 < 0 || y2 >= height) { continue; }
                    if (!bmap2[y2*width+x2]) { continue; }
                    matchable1[y1][x1] = 1;
                    matchable2[y2][x2] = 1;
                }
            }
        }
    }

    // Count the number of nodes on each side of the match.
    // Construct nodeID->pixel and pixel->nodeID maps.
    // Node IDs range from [0,n1) and [0,n2).
    int n1=0, n2=0;
    std::vector<Pixel> nodeToPix1;
    std::vector<Pixel> nodeToPix2;
    std::vector< std::vector< int > > pixToNode1 (height ,std::vector< int > (width, -1));
    std::vector< std::vector< int > > pixToNode2 (height ,std::vector< int > (width, -1));
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            Pixel pix (x,y);
            if (matchable1[y][x]) {
                pixToNode1[y][x] = n1;
                nodeToPix1.push_back(pix);
                n1++;
            }
            if (matchable2[y][x]) {
                pixToNode2[y][x] = n2;
                nodeToPix2.push_back(pix);
                n2++;
            }
        }
    }
    // Construct the list of edges between pixels within maxDist.
    std::vector<Edge> edges;
    for (int x1 = 0; x1 < width; ++x1) {
        for (int y1 = 0; y1 < height; ++y1) {
            if ( matchable1[y1][x1]==0) { continue; }
            for (int u = -r; u <= r; ++u) {
                for (int v = -r; v <= r; ++v) {
                    const double d2 = u*u + v*v;
                    if (d2 > maxDist*maxDist) { continue; }
                    const int x2 = x1 + u;
                    const int y2 = y1 + v;
                    if (x2 < 0 || x2 >= width) { continue; }
                    if (y2 < 0 || y2 >= height) { continue; }
                    if ( matchable2[y2][x2]==0 ) { continue; }
                    Edge e; 
                    e.i = pixToNode1[y1][x1];
                    e.j = pixToNode2[y2][x2];
                    e.w = sqrt(d2);
                    // printf("[%d,%d][%d,%d] %d %d | %d %d\n",matchable1[y1][x1],matchable2[y2][x2],e.i,e.j);
                    assert (e.i >= 0 && e.i < n1);
                    assert (e.j >= 0 && e.j < n2);
                    assert (e.w < outlierCost);
                    edges.push_back(e);
                }
            }
        }
    }

    // The cardinality of the match is n.
    const int n = n1 + n2;
    const int nmin = std::min(n1,n2);
    const int nmax = std::max(n1,n2);

    // Compute the degree of various outlier connections.
    const int d1 = std::max(0,std::min(degree,n1-1)); // from map1
    const int d2 = std::max(0,std::min(degree,n2-1)); // from map2
    const int d3 = std::min(degree,std::min(n1,n2)); // between outliers
    const int dmax = std::max(d1,std::max(d2,d3));

    assert (n1 == 0 || (d1 >= 0 && d1 < n1));
    assert (n2 == 0 || (d2 >= 0 && d2 < n2));
    assert (d3 >= 0 && d3 <= nmin);

    // Count the number of edges.
    int m = 0;
    m += edges.size();  // real connections
    m += d1 * n1;   // outlier connections
    m += d2 * n2;   // outlier connections
    m += d3 * nmax; // outlier-outlier connections
    m += n;         // high-cost perfect match overlay

    // If the graph is empty, then there's nothing to do.
    if (m == 0) {
        return 0;
    }

    // Weight of outlier connections.
    const int ow = (int) ceil (outlierCost * multiplier);

    // Scratch array for outlier edges.
    Array1D<int> outliers (dmax);

    // Construct the input graph for the assignment problem.
    Array2D<int> igraph (m,3);
    int count = 0;
    // real edges
    for (int a = 0; a < (int)edges.size(); ++a) {
        int i = edges[a].i;
        int j = edges[a].j;
        assert (i >= 0 && i < n1);
        assert (j >= 0 && j < n2);
        igraph(count,0) = i;
        igraph(count,1) = j;
        igraph(count,2) = (int) rint (edges[a].w * multiplier);
        count++;
    }
    // outliers edges for map1, exclude diagonal
    for (int i = 0; i < n1; ++i) {
        kOfN(d1,n1-1,outliers.data());
        for (int a = 0; a < d1; ++a) {
            int j = outliers(a);
            if (j >= i) { j++; }
            assert (i != j);
            assert (j >= 0 && j < n1);
            igraph(count,0) = i;
            igraph(count,1) = n2 + j;
            igraph(count,2) = ow;
            count++;
        }
    }
    // outliers edges for map2, exclude diagonal
    for (int j = 0; j < n2; ++j) {
        kOfN(d2,n2-1,outliers.data());
        for (int a = 0; a < d2; ++a) {
            int i = outliers(a);
            if (i >= j) { i++; }
            assert (i != j);
            assert (i >= 0 && i < n2);
            igraph(count,0) = n1 + i;
            igraph(count,1) = j;
            igraph(count,2) = ow;
            count++;
        }
    }
    // outlier-to-outlier edges
    for (int i = 0; i < nmax; ++i) {
        kOfN(d3,nmin,outliers.data());
        for (int a = 0; a < d3; ++a) {
            const int j = outliers(a);
            assert (j >= 0 && j < nmin);
            if (n1 < n2) {
                assert (i >= 0 && i < n2);
                assert (j >= 0 && j < n1);
                igraph(count,0) = n1 + i;
                igraph(count,1) = n2 + j;
            } else {
                assert (i >= 0 && i < n1);
                assert (j >= 0 && j < n2);
                igraph(count,0) = n1 + j;
                igraph(count,1) = n2 + i;
            }
            igraph(count,2) = ow;
            count++;
        }
    }
    // perfect match overlay (diagonal)
    for (int i = 0; i < n1; ++i) {
        igraph(count,0) = i;
        igraph(count,1) = n2 + i;
        igraph(count,2) = ow * multiplier;
        count++;
    }
    for (int i = 0; i < n2; ++i) {
        igraph(count,0) = n1 + i;
        igraph(count,1) = i;
        igraph(count,2) = ow * multiplier;
        count++;
    }
    assert (count == m);

    // Check all the edges, and set the values up for CSA.
    for (int i = 0; i < m; ++i) {
        assert(igraph(i,0) >= 0 && igraph(i,0) < n);
        assert(igraph(i,1) >= 0 && igraph(i,1) < n);
        igraph(i,0) += 1;
        igraph(i,1) += 1+n;
    }

    // Solve the assignment problem.
    CSA csa(2*n,m,igraph.data());
    assert(csa.edges()==n);

    Array2D<int> ograph (n,3);
    for (int i = 0; i < n; ++i) {
        int a,b,c;
        csa.edge(i,a,b,c);
        ograph(i,0)=a-1; ograph(i,1)=b-1-n; ograph(i,2)=c;
    }

    // Check the solution.
    // Count the number of high-cost edges from the perfect match
    // overlay that were used in the match.
    int overlayCount = 0;
    for (int a = 0; a < n; ++a) {
        const int i = ograph(a,0);
        const int j = ograph(a,1);
        const int c = ograph(a,2);
        assert (i >= 0 && i < n);
        assert (j >= 0 && j < n);
        assert (c >= 0);
        // edge from high-cost perfect match overlay
        if (c == ow * multiplier) { overlayCount++; }
        // skip outlier edges
        if (i >= n1) { continue; }
        if (j >= n2) { continue; }
        // for edges between real nodes, check the edge weight
        const Pixel pix1 = nodeToPix1[i];
        const Pixel pix2 = nodeToPix2[j];
        const int dx = pix1.x - pix2.x;
        const int dy = pix1.y - pix2.y;
        const int w = (int) rint (sqrt(dx*dx+dy*dy)*multiplier);
        assert (w == c);
    }

    // Print a warning if any of the edges from the perfect match overlay
    // were used.  This should happen rarely.  If it happens frequently,
    // then the outlier connectivity should be increased.
    if (overlayCount > 5) {
        fprintf (stderr, "%s:%d: WARNING: The match includes %d "
           "outlier(s) from the perfect match overlay.\n",
           __FILE__, __LINE__, overlayCount);
    }

    // Compute match arrays.
    for (int a = 0; a < n; ++a) {
        // node ids
        const int i = ograph(a,0);
        const int j = ograph(a,1);
        // skip outlier edges
        if (i >= n1) { continue; }
        if (j >= n2) { continue; }
        // map node ids to pixels
        const Pixel pix1 = nodeToPix1[i];
        const Pixel pix2 = nodeToPix2[j];
        // record edges
        match1[pix1.y][pix1.x] = pix2;
        match2[pix2.y][pix2.x] = pix1;
    }
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            if (bmap1[y*width+x]) {
                if (match1[y][x] != Pixel(-1,-1)) {
                    m1[y*width+x] = match1[y][x].y*width + match1[y][x].x;
                }
            }
            if (bmap2[y*width+x]) {
                if (match2[y][x] != Pixel(-1,-1)) {
                    m2[y*width+x] = match2[y][x].y*width + match2[y][x].x;
                }
            }
        }
    }

    // Compute the match cost.
    double cost = 0;
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            if (bmap1[y*width+x]) {
                if (match1[y][x] == Pixel(-1,-1)) {
                    cost += outlierCost;
                } else {
                    const int dx = x - match1[y][x].x;
                    const int dy = y - match1[y][x].y;
                    cost += 0.5 * sqrt (dx*dx + dy*dy);
                }
            }
            if (bmap2[y*width+x]) {
                if (match2[y][x] == Pixel(-1,-1)) {
                    cost += outlierCost;
                } else {
                    const int dx = x - match2[y][x].x;
                    const int dy = y - match2[y][x].y;
                    cost += 0.5 * sqrt (dx*dx + dy*dy);
                }
            }
        }
    }    

    // Return the match cost.
    return cost;
}


double matchEdgeMaps3D(const double* bmap1, const int height, const int width, const int depth, 
                     const double* bmap2, const int height2, const int width2, const int depth2,
                     double* m1, int m1h, int m1w, int m1z,
                     double* m2, int m2h, int m2w, int m2z,
                     double maxDist, double outlierCost, int degree)
{

    assert (degree > 0);
    assert (multiplier > 0);

    // Check arguments.
    assert (height == height2);
    assert (width == width2);
    assert (depth == depth2);

    assert (maxDist >= 0);
    assert (outlierCost > maxDist);

    // Initialize match[12] arrays to (-1,-1,-1).
    Array3D<Voxel> match1 (width,height,depth);
    Array3D<Voxel> match2 (width,height,depth);
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            for (int z = 0; z < depth; ++z) {
                match1(x,y,z) = Voxel(-1,-1,-1);
                match2(x,y,z) = Voxel(-1,-1,-1);
            }
        }
    }

    // Radius of search window.
    const int r = (int) ceil (maxDist);
    printf("max dist: %f \n", maxDist);

    // Figure out which nodes are matchable, i.e. within maxDist
    // of another node.
    // printf("%s\n","Figure out which nodes are matchable.");
    Array3D<bool> matchable1 (width,height,depth);
    Array3D<bool> matchable2 (width,height,depth);
    matchable1.init(false);
    matchable2.init(false);
    for (int y1 = 0; y1 < height; ++y1) {
        for (int x1 = 0; x1 < width; ++x1) {
            for (int z1 = 0; z1 < depth; ++z1) {
                if (!bmap1[(y1*width+x1)*depth+z1]) { continue; }
                for (int v = -r; v <= r; ++v) {
                    for (int u = -r; u <= r; ++u) {
                        for (int w = -r; w <= r; ++w) {
                            const double d2 = u*u + v*v + w*w;
                            if (d2 > maxDist*maxDist) { continue; }
                            const int x2 = x1 + u;
                            const int y2 = y1 + v;
                            const int z2 = z1 + w;
                            if (x2 < 0 || x2 >= width) { continue; }
                            if (y2 < 0 || y2 >= height) { continue; }
                            if (z2 < 0 || z2 >= depth) { continue; }
                            if (!bmap2[(y1*width+x1)*depth+z1]) { continue; }
                            matchable1(x1,y1,z1) = true;
                            matchable2(x2,y2,z2) = true;
                        }
                    }
                }
            }
        }
    }

    // Count the number of nodes on each side of the match.
    // Construct nodeID->pixel and pixel->nodeID maps.
    // Node IDs range from [0,n1) and [0,n2).
    int n1=0, n2=0;
    std::vector<Voxel> nodeToPix1;
    std::vector<Voxel> nodeToPix2;
    // printf("%s\n","Count the number of nodes on each side of the match.");
    Array3D<int> pixToNode1 (width,height,depth);
    Array3D<int> pixToNode2 (width,height,depth);
    
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            for (int z = 0; z < depth; ++z) {
                pixToNode1(x,y,z) = -1;
                pixToNode2(x,y,z) = -1;
                Voxel pix (x,y,z);
                if (matchable1(x,y,z)) {
                    pixToNode1(x,y,z) = n1;
                    nodeToPix1.push_back(pix);
                    n1++;
                }
                if (matchable2(x,y,z)) {
                    pixToNode2(x,y,z) = n2;
                    nodeToPix2.push_back(pix);
                    n2++;
                }
            }
        }
    }
    // printf("done %s\n","Count the number of nodes on each side of the match.");

    // Construct the list of edges between pixels within maxDist.
    // printf("%s\n","Construct the list of edges between pixels within maxDist.");
    std::vector<Edge> edges;
    for (int x1 = 0; x1 < width; ++x1) {
        for (int y1 = 0; y1 < height; ++y1) {
            for (int z1 = 0; z1 < depth; ++z1) {
                if (!matchable1(x1,y1,z1)) { continue; }
                for (int u = -r; u <= r; ++u) {
                    for (int v = -r; v <= r; ++v) {
                        for (int w = -r; w <= r; ++w) {
                            const double d2 = u*u + v*v + w*w;
                            if (d2 > maxDist*maxDist) { continue; }
                            const int x2 = x1 + u;
                            const int y2 = y1 + v;
                            const int z2 = z1 + w;
                            if (x2 < 0 || x2 >= width) { continue; }
                            if (y2 < 0 || y2 >= height) { continue; }
                            if (z2 < 0 || z2 >= depth) { continue; }
                            if (!matchable2(x2,y2,z2)) { continue; }
                            Edge e; 
                            e.i = pixToNode1(x1,y1,z1);
                            e.j = pixToNode2(x2,y2,z2);
                            e.w = sqrt(d2);
                            assert (e.i >= 0 && e.i < n1);
                            assert (e.j >= 0 && e.j < n2);
                            assert (e.w < outlierCost);
                            edges.push_back(e);
                        }
                    }
                }
            }
        }
    }
   // printf("done %s\n","Construct the list of edges between pixels within maxDist.");

    // The cardinality of the match is n.
    const int n = n1 + n2;
    const int nmin = std::min(n1,n2);
    const int nmax = std::max(n1,n2);

    // Compute the degree of various outlier connections.
    const int d1 = std::max(0,std::min(degree,n1-1)); // from map1
    const int d2 = std::max(0,std::min(degree,n2-1)); // from map2
    const int d3 = std::min(degree,std::min(n1,n2)); // between outliers
    const int dmax = std::max(d1,std::max(d2,d3));

    assert (n1 == 0 || (d1 >= 0 && d1 < n1));
    assert (n2 == 0 || (d2 >= 0 && d2 < n2));
    assert (d3 >= 0 && d3 <= nmin);

    // Count the number of edges.
    int m = 0;
    m += edges.size();  // real connections
    m += d1 * n1;   // outlier connections
    m += d2 * n2;   // outlier connections
    m += d3 * nmax; // outlier-outlier connections
    m += n;         // high-cost perfect match overlay

    // If the graph is empty, then there's nothing to do.
    if (m == 0) {
        return 0;
    }

    // Weight of outlier connections.
    const int ow = (int) ceil (outlierCost * multiplier);

    // Scratch array for outlier edges.
    Array1D<int> outliers (dmax);

    // Construct the input graph for the assignment problem.
    // igraph is (node1, node2, weight)
    Array2D<int> igraph (m,3);
    int count = 0;
    // real edges
    for (int a = 0; a < (int)edges.size(); ++a) {
        int i = edges[a].i;
        int j = edges[a].j;
        assert (i >= 0 && i < n1);
        assert (j >= 0 && j < n2);
        igraph(count,0) = i;
        igraph(count,1) = j;
        igraph(count,2) = (int) rint (edges[a].w * multiplier);
        count++;
    }

    // outliers edges for map1, exclude diagonal
    for (int i = 0; i < n1; ++i) {
        kOfN(d1,n1-1,outliers.data());
        for (int a = 0; a < d1; ++a) {
            int j = outliers(a);
            if (j >= i) { j++; }
            assert (i != j);
            assert (j >= 0 && j < n1);
            igraph(count,0) = i;
            igraph(count,1) = n2 + j;
            igraph(count,2) = ow;
            count++;
        }
    }
    // outliers edges for map2, exclude diagonal
    for (int j = 0; j < n2; ++j) {
        kOfN(d2,n2-1,outliers.data());
        for (int a = 0; a < d2; ++a) {
            int i = outliers(a);
            if (i >= j) { i++; }
            assert (i != j);
            assert (i >= 0 && i < n2);
            igraph(count,0) = n1 + i;
            igraph(count,1) = j;
            igraph(count,2) = ow;
            count++;
        }
    }
    // outlier-to-outlier edges
    for (int i = 0; i < nmax; ++i) {
        kOfN(d3,nmin,outliers.data());
        for (int a = 0; a < d3; ++a) {
            const int j = outliers(a);
            assert (j >= 0 && j < nmin);
            if (n1 < n2) {
                assert (i >= 0 && i < n2);
                assert (j >= 0 && j < n1);
                igraph(count,0) = n1 + i;
                igraph(count,1) = n2 + j;
            } else {
                assert (i >= 0 && i < n1);
                assert (j >= 0 && j < n2);
                igraph(count,0) = n1 + j;
                igraph(count,1) = n2 + i;
            }
            igraph(count,2) = ow;
            count++;
        }
    }
    // perfect match overlay (diagonal)
    for (int i = 0; i < n1; ++i) {
        igraph(count,0) = i;
        igraph(count,1) = n2 + i;
        igraph(count,2) = ow * multiplier;
        count++;
    }
    for (int i = 0; i < n2; ++i) {
        igraph(count,0) = n1 + i;
        igraph(count,1) = i;
        igraph(count,2) = ow * multiplier;
        count++;
    }
    assert (count == m);

    // Check all the edges, and set the values up for CSA.
    for (int i = 0; i < m; ++i) {
        assert(igraph(i,0) >= 0 && igraph(i,0) < n);
        assert(igraph(i,1) >= 0 && igraph(i,1) < n);
        igraph(i,0) += 1;
        igraph(i,1) += 1+n;
    }

    // Solve the assignment problem.
    CSA csa(2*n,m,igraph.data());
    assert(csa.edges()==n);

    Array2D<int> ograph (n,3);
    for (int i = 0; i < n; ++i) {
        int a,b,c;
        csa.edge(i,a,b,c);
        ograph(i,0)=a-1; ograph(i,1)=b-1-n; ograph(i,2)=c;
    }

    // Check the solution.
    // Count the number of high-cost edges from the perfect match
    // overlay that were used in the match.
    int overlayCount = 0;
    for (int a = 0; a < n; ++a) {
        const int i = ograph(a,0);
        const int j = ograph(a,1);
        const int c = ograph(a,2);
        assert (i >= 0 && i < n);
        assert (j >= 0 && j < n);
        assert (c >= 0);
        // edge from high-cost perfect match overlay
        if (c == ow * multiplier) { overlayCount++; }
        // skip outlier edges
        if (i >= n1) { continue; }
        if (j >= n2) { continue; }
        // for edges between real nodes, check the edge weight
        const Voxel pix1 = nodeToPix1[i];
        const Voxel pix2 = nodeToPix2[j];
        const int dx = pix1.x - pix2.x;
        const int dy = pix1.y - pix2.y;
        const int dz = pix1.z - pix2.z;
        const int w = (int) rint (sqrt(dx*dx+dy*dy+dz*dz)*multiplier);
        assert (w == c);
    }

    // Print a warning if any of the edges from the perfect match overlay
    // were used.  This should happen rarely.  If it happens frequently,
    // then the outlier connectivity should be increased.
    if (overlayCount > 5) {
        fprintf (stderr, "%s:%d: WARNING: The match includes %d "
         "outlier(s) from the perfect match overlay.\n",
         __FILE__, __LINE__, overlayCount);
    }

    // Compute match arrays.
    for (int a = 0; a < n; ++a) {
        // node ids
        const int i = ograph(a,0);
        const int j = ograph(a,1);
        // skip outlier edges
        if (i >= n1) { continue; }
        if (j >= n2) { continue; }
        // map node ids to pixels
        const Voxel pix1 = nodeToPix1[i];
        const Voxel pix2 = nodeToPix2[j];
        // record edges
        match1(pix1.x,pix1.y,pix1.z) = pix2;
        match2(pix2.x,pix2.y,pix2.z) = pix1;
    }
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            for (int z = 0; z < depth; ++z) {
                if (bmap1[(y*width+x)*depth+z]) {
                    if (match1(x,y,z) != Voxel(-1,-1,-1)) {
                        m1[(y*width+x)*depth+z] = 
                        match1(x,y,z).y*width*depth + 
                        match1(x,y,z).x*depth +
                        match1(x,y,z).z +
                        1;
                    }
                }
                if (bmap2[(y*width+x)*depth+z]) {
                    if (match2(x,y,z) != Voxel(-1,-1,-1)) {
                        m2[(y*width+x)*depth+z] = 
                        match2(x,y,z).y*width*depth + 
                        match2(x,y,z).x*depth +
                        match2(x,y,z).z +
                        1;
                    }
                }
            }
        }
    }

    // Compute the match cost.
    double cost = 0;
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            for (int z = 0; z < depth; ++z) {
                if (bmap1[(y*width+x)*depth+z]) {
                    if (match1(x,y,z) == Voxel(-1,-1,-1)) {
                        cost += outlierCost;
                    } else {
                        const int dx = x - match1(x,y,z).x;
                        const int dy = y - match1(x,y,z).y;
                        const int dz = z - match1(x,y,z).z;
                        cost += 0.5 * sqrt (dx*dx + dy*dy + dz*dz);
                    }
                }
                if (bmap2[(y*width+x)*depth+z]) {
                    if (match2(x,y,z) == Voxel(-1,-1,-1)) {
                        cost += outlierCost;
                    } else {
                        const int dx = x - match2(x,y,z).x;
                        const int dy = y - match2(x,y,z).y;
                        const int dz = z - match2(x,y,z).z;
                        cost += 0.5 * sqrt (dx*dx + dy*dy + dz*dz);
                    }
                }
            }
        }
    }

    // Return the match cost.
    return cost;
}


// double
// matchEdgeMaps (
//     Matrix& bmap1, Matrix& bmap2,
//     double maxDist, double outlierCost,
//     Matrix& m1, Matrix& m2)
// {
//     m1 = Matrix(bmap1.nrows(),bmap1.ncols());
//     m2 = Matrix(bmap2.nrows(),bmap2.ncols());
//     return matchEdgeMaps2D(bmap1.data(), bmap1.nrows(), bmap1.ncols(), 
//       bmap2.data(), bmap2.nrows(), bmap2.ncols(), 
//       m1.data(), m1.nrows(), m1.ncols(), 
//       m2.data(), m2.nrows(), m2.ncols(),
//       maxDist, outlierCost);

// }
