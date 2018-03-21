#include "poly_tri.h"
// TODO: 
//  - test delaunay
//  - test hole
//  - find problem in triangle creation

void test_0()
{
    Vecs pts
    {
        Vec{-1., 0.},
        Vec{1.,  0.},
        Vec{0.,  0.5},
        Vec{0., -0.5}
    };
    Boundaries b {Boundary{0, 1}};
    PolyTri p(pts, b); 
    std::cout << "\nnum of tris: "  << p.tris.size() << std::endl;
    std::cout << "size of point2tris: " << p.pnt2tris.size() << std::endl;
    for (auto i: p.get_tris())
    {
        for (auto j: i)
            std::cout << j;
        std::cout << std::endl;
    }
    for (int i=0; i<4; i++)
        std::cout << p.pnt2tris[i].size() << std::endl;
}

void test_1()
{
    Vecs pts
    {
        Vec{-1., 0.},
        Vec{1.,  0.},
        Vec{0.5,  0.5},
        Vec{0.5, -0.5},
        Vec{-0.5,  0.5},
        Vec{-0.5, -0.5}
    };
    Boundaries b {Boundary{0, 1}};
    PolyTri p(pts, b); 
    std::cout << "\nnum of tris: "  << p.tris.size() << std::endl;
    std::cout << "size of point2tris: " << p.pnt2tris.size() << std::endl;
    for (auto i: p.get_tris())
    {
        for (auto j: i)
            std::cout << j;
        std::cout << std::endl;
    }
    for (int i=0; i<4; i++)
        std::cout << p.pnt2tris[i].size() << std::endl;
}



void test_2()
{
    Vecs pts
    {
        Vec{0, 0},
        Vec{0.2, 0.1},
        Vec{0.5, 0.1},
        Vec{0.8, 0.1},
        Vec{1.,  0. }
    };
    Boundaries b {Boundary{1, 2, 3, 4, 0, 1}};
    PolyTri p(pts, b); 
    std::cout << "\nnum of tris: "  << p.tris.size() << std::endl;
    std::cout << "size of point2tris: " << p.pnt2tris.size() << std::endl;
    for (auto i: p.get_tris())
    {
        for (auto j: i)
            std::cout << j;
        std::cout << std::endl;
    }
    for (int i=0; i<4; i++)
        std::cout << p.pnt2tris[i].size() << std::endl;
}


void test_3()
{
    Vecs pts{
        Vec{ 0.        ,  0.  },
        Vec{ 0.11111111, -0.1 },
        Vec{ 0.22222222, -0.1 },
        Vec{ 0.33333333, -0.1 },
        Vec{ 0.44444444, -0.1 },
        Vec{ 0.55555556, -0.1 },
        Vec{ 0.66666667, -0.1 },
        Vec{ 0.77777778, -0.1 },
        Vec{ 0.88888889, -0.1 },
        Vec{ 1.        ,  0.  },
        Vec{ 0.11111111,  0.1 },
        Vec{ 0.22222222,  0.1 },
        Vec{ 0.33333333,  0.1 },
        Vec{ 0.44444444,  0.1 },
        Vec{ 0.55555556,  0.1 },
        Vec{ 0.66666667,  0.1 },
        Vec{ 0.77777778,  0.1 },
        Vec{ 0.88888889,  0.1 },
        Vec{-1.        ,  0.  },
        Vec{ 2.        ,  0.  }
    };
    Boundaries b {Boundary{0, 9}};
    PolyTri p(pts, b);
    std::cout << "\nnum of tris: "  << p.tris.size() << std::endl;
    std::cout << "size of point2tris: " << p.pnt2tris.size() << std::endl;
}

void test_4()
{
    Vecs pts{
        Vec{ 0.        ,  0.  },
        Vec{ 1.        ,  0.  },
        Vec{ 0.11111111,  0.1 },
        Vec{ 0.22222222,  0.1 },
        Vec{ 0.33333333,  0.1 },
        Vec{ 0.44444444,  0.1 },
        Vec{ 0.55555556,  0.1 },
        Vec{ 0.66666667,  0.1 },
        Vec{ 0.77777778,  0.1 },
        Vec{ 0.88888889,  0.1 }
    };
    Boundaries b {Boundary{0, 1}};
    PolyTri p(pts, b);
    std::cout << "\nnum of tris: "  << p.tris.size() << std::endl;
    for (auto i: p.get_tris())
    {
        for (auto j: i)
            std::cout << j;
        std::cout << std::endl;
    }
    for (int i=0; i<4; i++)
        std::cout << p.pnt2tris[i].size() << std::endl;
}


void test_5()
{
            Vecs pts{
                Vec{ 5.00000000e-01,  0.00000000e+00},
                Vec{-2.50000000e-01,  4.33012702e-01},
                Vec{-2.50000000e-01, -4.33012702e-01},
                Vec{ 5.00000000e-01, -1.22464680e-16},
                Vec{ 2.00000000e+00,  0.00000000e+00},
                Vec{-1.00000000e+00,  8.66025404e-01},
                Vec{-1.00000000e+00, -8.66025404e-01},
                Vec{ 2.00000000e+00, -2.44929360e-16}
            };
            Boundaries b {Boundary{0, 1, 2, 3}};
            PolyTri p(pts, b);
}


int main()
{
    test_5();
    return 0;
}
