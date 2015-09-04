#ifndef MESH_H
#define MESH_H

template <typename D>
class Mesh
{
public:
   Mesh();
   ~Mesh();

   void create_uniform(int num_intervals);
   void create_refined(int num_intervals);

   void interval_endpoints(int interval_idx, D &a, D &b);

   // as degrees of freedom: idx = 0 -> first internal point (not 0) is returned
   D gridpoint(int internal_gridpoint_idx);

   int num_intervals;
   int N;

public: //private:
   void allocate(int num_intervals);

private:
   D *gridpoints;
};


template <typename D>
void Mesh<D>::interval_endpoints(int interval_idx, D &a, D &b)
{
   assert((interval_idx >= 0) && (interval_idx < this->num_intervals));
   if(interval_idx == 0)
   {
      a = 0;
   }
   else
   {
      a = gridpoints[interval_idx - 1];
   }

   if(interval_idx == this->num_intervals - 1)
   {
      b = 1;
   }
   else
   {
      b = gridpoints[interval_idx];
   }

}

template <typename D>
void Mesh<D>::create_uniform(int num_intervals)
{
   allocate(num_intervals);

   D h = 1./num_intervals;

   for(int i = 0; i < N; i++)
      gridpoints[i] = h*(i+1);
}

template <typename D>
void Mesh<D>::create_refined(int num_intervals)
{
   allocate(num_intervals);

   D q = pow(N, 1./(N-1.));
   D a = (1-q)/(1-pow(q,N));
   D sum = 0.;
   for(int i = 0; i < N; i++)
   {
      D int_len = a * pow(q, i);
      sum += int_len;
      gridpoints[i] = sum;
   }
}


//void distribute(int N)
//{
//   const double q = pow(N, 1./(N-1.));
//   double a = (1-q)/(1-pow(q,N));
//   double sum = 0;
//   for(int i = 0; i < N; i++)
//   {
//      double x = a * pow(q, i);
//      sum += x;
//      cout << x << endl;// ";   " << 1 - x << endl;
//   }
//   cout << "sum " << sum << ", first " << a << ",last " << a*pow(q, N-1) << ", first/last ratio " << pow(q, N-1) << endl;
//}


template <typename D>
Mesh<D>::Mesh()
{
   gridpoints = nullptr;
}

template <typename D>
Mesh<D>::~Mesh()
{
   if(gridpoints == nullptr)
      delete[] gridpoints;
}

template <typename D>
void Mesh<D>::allocate(int num_intervals)
{
   if(gridpoints == nullptr)
      delete[] gridpoints;

   this->num_intervals = num_intervals;
   N = num_intervals - 1;

   gridpoints = new D[N];
}

template <typename D_in, typename D_out>
void convert_mesh(const Mesh<D_in> &mesh_in, Mesh<D_out> &mesh_out)
{
   mesh_out.N = mesh_in.N;
   mesh_out.num_intervals = mesh_in.num_intervals;
   mesh_out.allocate(mesh_in.num_intervals);
}


#endif // MESH_H
