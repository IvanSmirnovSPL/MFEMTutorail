#include "mfem.hpp"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

double f1(const Vector &x)
{

      return cos(2 * M_PI * x(0)) * cos(2 * M_PI * x(1));
}

double f0(const Vector &x)
{

      return 1;
}

int main(int argc, char *argv[])
{
    // 1. Parse command-line options.
    int ref_levels = 0, rightPart = 0;
    OptionsParser args(argc, argv);
    args.AddOption(&ref_levels, "-refine_levels", "--refine_levels",
                  "The number of mesh refinements.");
    args.AddOption(&rightPart, "-rightPart", "--rightPart",
     "Right part: 0 -> f = 1, 1 -> f = cos(2 * Pi * x) * cos( 2 * Pi * y)");
                  
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }

   double(*f)(const Vector &x);
   if (rightPart == 0)
   {
       f = &f0;
   }
   else
   {
       f = &f1;
   }

    const char *mesh_file = "./star.mesh";
    int order = 1;
   // bool static_cond = false;
//    bool pa = false;
//    bool fa = false;
//    const char *device_config = "cpu";
//    bool visualization = true;
//    bool algebraic_ceed = false;

    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    
    int dim = mesh->Dimension();
    for(int i = 0; i < ref_levels; i += 1)
        mesh->UniformRefinement();

    std::cout<<"number of elements: "<<mesh->GetNE()<<std::endl;
//    {
//        int ref_levels =
//                (int)floor(log(50000./mesh->GetNE())/log(2.)/dim);
//        for (int l = 0; l < ref_levels; l++)
//        {
//            mesh->UniformRefinement();
//        }
//    }

    FiniteElementCollection *fec = new H1_FECollection(order, dim);
    FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);

    Array<int> ess_tdof_list;
    Array<int> ess_bdr;
    ess_bdr.SetSize(mesh->bdr_attributes.Max());
    ess_bdr = 1;
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    LinearForm b(fespace);
    ConstantCoefficient one(1.0);
    FunctionCoefficient rhs_coef(f);
    //b.AddDomainIntegrator(new DomainLFIntegrator(one));
    b.AddDomainIntegrator(new DomainLFIntegrator(rhs_coef));
    b.Assemble();

    GridFunction x(fespace);
    x = 0.0;

    BilinearForm a(fespace);
    a.AddDomainIntegrator(new DiffusionIntegrator(one));
    a.Assemble();

    OperatorPtr A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

//    GSSmoother M((SparseMatrix&)(*A));
//    PCG(*A, M, B, X, 1, 200, 1e-12, 0.0);

    UMFPackSolver umf_solver;
    umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
    umf_solver.SetOperator(*A);
    umf_solver.Mult(B, X);

    a.RecoverFEMSolution(X, b, x);
    string name = "_"+ std::to_string(ref_levels) +"_" +std::to_string(rightPart);

    ofstream mesh_ofs("refined"+ name +".mesh");
    mesh_ofs.precision(8);
    mesh->Print(mesh_ofs);
    ofstream sol_ofs("sol"+ name +".gf");
    sol_ofs.precision(8);
    x.Save(sol_ofs);

    ofstream vtk_ofs("out_sol"+ name +".vtk");
    vtk_ofs.precision(8);
    int ref = 0;
    mesh->PrintVTK(vtk_ofs, ref);
    x.SaveVTK(vtk_ofs, "u", ref);
    vtk_ofs.close();


    return 0;
}
