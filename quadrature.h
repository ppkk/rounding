#ifndef QUADRATURE_H
#define QUADRATURE_H


template <typename R>
class Quadrature
{
public:
   Quadrature(int np);
   ~Quadrature();
   int n_points;
   R* weights;
   R* points;
};

template <typename R>
Quadrature<R>::Quadrature(int np)
{
   n_points = np;
   if(np == 3)
   {
      weights = new R[np] {5./9., 8./9., 5./9.};
      points = new R[np] {-sqrt(3./5.), 0, sqrt(3./5.)};
   }
   else
      assert(0);
}

template <typename R>
Quadrature<R>::~Quadrature()
{
   delete[] weights;
   delete[] points;
}


#endif // QUADRATURE_H
