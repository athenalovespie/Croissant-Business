#include "energy.h"
#include "tests.h"
#include "houdini_export.h"
#include <cmath>
#include <iostream>
#include <tuple>
#include <vector>

std::tuple<Eigen::VectorXd, Eigen::MatrixXd, Eigen::MatrixXd> get_svd(const Matrix32& F){
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeThinU | Eigen::ComputeThinV);
    return {
        svd.singularValues(), svd.matrixU(), svd.matrixV()
    };
}

double psi(const Matrix32 &F, double lambda, double mew){

    Eigen::JacobiSVD<Matrix32> svd(F);
    Eigen::VectorXd sigmas = svd.singularValues();
    double sigma_1 = sigmas(0);
    double sigma_2 = sigmas(1);
    double E = (sigma_1*sigma_1 + sigma_2*sigma_2)-2*(sigma_1+sigma_2) + 2;
    double P = ((sigma_1*sigma_2) - 1)*((sigma_1*sigma_2) - 1);
    return mew*E + .5*lambda*P;
}

Matrix32 D_psi(const Matrix32& F, double lambda, double mew){


    auto [sigmas, U, V] = get_svd(F);
    double sigma_1 = sigmas(0);
    double sigma_2 = sigmas(1);


    Eigen::MatrixXd R = U*V.transpose();
    Eigen::MatrixXd C_inv = (F.transpose()*F).inverse();
    double detFTF = sigma_1*sigma_1*sigma_2*sigma_2;
    return mew*(2*F -2*R) + F*C_inv*lambda*(detFTF - std::sqrt(detFTF));
}

Matrix32 DD_psi(const Matrix32& F, const Matrix32& delta_F, double lambda, double mew){
    Matrix33 I = Matrix33::Identity();
    auto [sigmas, U, V] = get_svd(F);
    double sigma_1 = sigmas(0);
    double sigma_2 = sigmas(1);

    Eigen::MatrixXd R = U*V.transpose();
    Eigen::MatrixXd C_inv = (F.transpose()*F).inverse();
    Matrix22 Sigma;
    Sigma << sigma_1, 0.0, 0.0, sigma_2;

    Matrix22 S = V * Sigma * V.transpose();

    double J_squared = sigma_1*sigma_1*sigma_2*sigma_2;
    double J = sigma_1*sigma_2;

    return
    2.0*mew*(
        delta_F
        - (1.0/(sigma_1 + sigma_2))
            * R
            * (R.transpose()*delta_F - delta_F.transpose()*R)
        - (I - R*R.transpose())
            * delta_F
            * S.inverse()
    )
    + lambda*(
        J*(2.0*J - 1.0)
            * ((F*C_inv).array()*delta_F.array()).sum()
            * F*C_inv

        + (J_squared - J)
            * delta_F*C_inv

        - (J_squared - J)
            * F*C_inv
            * (delta_F.transpose()*F + F.transpose()*delta_F)
            * C_inv
    );

}


double area(int t, Eigen::VectorXi& mesh, Eigen::VectorXd& reference_components, Eigen::VectorXd& current_components){
    int V1 = mesh(3*t);
    int V2 = mesh(3*t+1);
    int V3 = mesh(3*t+2);

    double x11 = reference_components(3*V1);
    double x12 = reference_components(3*V1+1);

    double x21 = reference_components(3*V2);
    double x22 = reference_components(3*V2+1);

    double x31 = reference_components(3*V3);
    double x32 = reference_components(3*V3+1);

    return .5*std::abs(x11*(x22 - x32) + x21*(x32 - x12) +x31*(x12 -x22));
}

Eigen::VectorXd get_vertexmass(Eigen::VectorXi& mesh, Eigen::VectorXd& reference_components, Eigen::VectorXd& current_components){
    double p = 1000;
    double h = 0.001;
    
    Eigen::VectorXd masses = Eigen::VectorXd::Zero(reference_components.size() / 3);
    for (int t = 0; t < mesh.size()/3; t++){
        int V1 = mesh(3*t);
        int V2 = mesh(3*t+1);
        int V3 = mesh(3*t+2);

        double mass = p*area(t, mesh, reference_components, current_components)*h;
        masses(V1) += mass/3;
        masses(V2) += mass/3;
        masses(V3) += mass/3;
    }
    return masses;
}


Matrix32 DS(int t, Eigen::VectorXi& mesh, Eigen::VectorXd& reference_components, Eigen::VectorXd& current_components){

    int vertex1 = mesh(3*t);
    int vertex2 = mesh(3*t+1);
    int vertex3 = mesh(3*t+2);

    //vertex 1 of triangle components
    double x11 = current_components(3*vertex1);
    double x12 = current_components(3*vertex1+1);
    double x13 = current_components(3*vertex1+2);

    double x21 = current_components(3*vertex2);
    double x22 = current_components(3*vertex2+1);
    double x23 = current_components(3*vertex2+2);

    double x31 = current_components(3*vertex3);
    double x32 = current_components(3*vertex3+1);
    double x33 = current_components(3*vertex3+2);

    // for(int i = 0; i < 9*N; i++){
    // Eigen::VectorXi& vertices = mesh(3*t + (i-1));
    // Eigen::VectorXd& x = current_components(3*vertices(i) + (i-1)); //x(i)
    // }
    //vertex 2 of triangle components

    Matrix32 DS;
    DS << (x11 - x31), (x21 - x31),
          (x12 - x32), (x22 - x32),
          (x13 - x33), (x23 - x33);


    // double fraction = 1/((x11-x31)*(x22-x32) - (x21 - x31)(x12 - x32));
    return DS;
}

Matrix22 Dm_inv(int t, Eigen::VectorXi& mesh, Eigen::VectorXd& reference_components, Eigen::VectorXd& current_components){
    int vertex1 = mesh(3*t);
    int vertex2 = mesh(3*t+1);
    int vertex3 = mesh(3*t+2);

    //vertex 1 of triangle components
    double X11 = reference_components(3*vertex1);
    double X12 = reference_components(3*vertex1+1);
    double X13 = reference_components(3*vertex1+2);

    double X21 = reference_components(3*vertex2);
    double X22 = reference_components(3*vertex2+1);
    double X23 = reference_components(3*vertex2+2);

    double X31 = reference_components(3*vertex3);
    double X32 = reference_components(3*vertex3+1);
    double X33 = reference_components(3*vertex3+2);

    Matrix22 Dm_inv;
    Dm_inv << (X22 -X32), (X31 - X21),
              (X32 - X12), (X11 - X31);

    Dm_inv*= 1/((X11 - X31)*(X22 - X32) - (X21 - X31)*(X12 - X32));

    return Dm_inv;
}

Matrix32 F(int t, Eigen::VectorXi& mesh, Eigen::VectorXd& reference_components, Eigen::VectorXd& current_components){
    Matrix32 F = DS(t,mesh,reference_components,current_components)*Dm_inv(t,mesh,reference_components,current_components);
    return F;
}

double PE(int N, Eigen::VectorXi& mesh, Eigen::VectorXd& reference_components, Eigen::VectorXd& current_components, double mew, double lambda){
    double PE = 0;
    for (int t = 0; t < N; ++t){
        PE += area(t,mesh,reference_components,current_components)*psi(F(t, mesh, reference_components, current_components), lambda, mew);
    }
    return PE;
}

Eigen::VectorXd gradient(int t, int N, Eigen::VectorXi& mesh, Eigen::VectorXd& reference_components, Eigen::VectorXd& current_components, double mew, double lambda){
    Eigen::VectorXd f = Eigen::VectorXd::Zero(current_components.size());


    double DFX;

    for(int t = 0; t < N; t++){
        Matrix32 F1 = F(t, mesh, reference_components, current_components);
        double area1 = area(t, mesh, reference_components, current_components);

        int vertex1 = mesh(3*t);
        int vertex2 = mesh(3*t+1);
        int vertex3 = mesh(3*t+2);

        double X11 = reference_components(3*vertex1);
        double X12 = reference_components(3*vertex1+1);
        double X13 = reference_components(3*vertex1+2);

        double X21 = reference_components(3*vertex2);
        double X22 = reference_components(3*vertex2+1);
        double X23 = reference_components(3*vertex2+2);

        double X31 = reference_components(3*vertex3);
        double X32 = reference_components(3*vertex3+1);
        double X33 = reference_components(3*vertex3+2);

        double c = 1/((X11 - X31)*(X22 - X32) - (X21 - X31)*(X12 - X32));
            for(int i = 0; i < 3; i++){
                Matrix32 P = D_psi(F1, lambda, mew);
                for(int a = 0; a < 3; a++){
                    for(int g = 0; g < 2; g ++){
                        if (g == 0 && i == 0){
                            DFX = c*(X22 -X32);
                        }
                        else if (g == 1 && i==0){
                            DFX = c*(X31-X21);
                        }
                        else if (g == 0 && i==1){
                            DFX = c*(X32 - X12);
                        }
                        else if (g == 1 && i==1){
                            DFX = c*(X11 - X31);
                        }
                        else if (g == 0 && i==2){
                            DFX = c*(X12-X22);
                        }
                        else if (g == 1 && i==2){
                            DFX = c*(X21 - X11);
                        }

                        f[3*mesh[3*t +i] + a] += P(a, g)*DFX*area1;


                }
            }
        }

    }
    return f;
}

Eigen::VectorXd grav_gradient(int t, int N, Eigen::VectorXi& mesh, Eigen::VectorXd& reference_components, Eigen::VectorXd& current_components, double mew, double lambda){
    Eigen::VectorXd masses = get_vertexmass(mesh, reference_components, current_components);
    Eigen::VectorXd G = Eigen::VectorXd::Zero(current_components.size());;
        
    for (int i = 0; i < masses.size(); i++){   
        G(3*i+1)= masses(i)*(-9.8);
    }
    
    return gradient(0, N, mesh, reference_components, current_components, mew, lambda) - G;
}

Eigen::VectorXd gradient_descent(int N,
    Eigen::VectorXi& mesh,
    Eigen::VectorXd& reference_components,
    Eigen::VectorXd& current_components,
    double mew,
    double lambda,
    const std::vector<int>& pinned_vertices){
        
    Eigen::VectorXd x = current_components;

    Eigen::VectorXd masses =
        get_vertexmass(mesh, reference_components, x);

    // Total potential energy: elasticity + gravity.
    auto total_energy = [&](Eigen::VectorXd& positions) {
        double energy =
            PE(N, mesh, reference_components, positions, mew, lambda);

        for (Eigen::Index i = 0; i < masses.size(); ++i) {
            energy += 9.8 * masses(i) * positions(3*i + 1);
        }

        return energy;
    };

    const int max_iterations = 1000;
    const double tolerance = 1e-6;

    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        Eigen::VectorXd g =
            grav_gradient(
                0, N, mesh, reference_components, x, mew, lambda
            );

        // Pinned vertices must not move in any coordinate.
        for (int vertex : pinned_vertices) {
            g.segment<3>(3*vertex).setZero();
        }

        if (!g.allFinite()) {
            throw std::runtime_error("Gradient contains invalid values.");
        }

        // Stop when the gradient on the free vertices is small.
        if (g.norm() < tolerance) {
            return x;
        }

        double energy_before = total_energy(x);
        double alpha = 0.1;
        bool accepted = false;

        for (int attempt = 0; attempt < 40; ++attempt) {
            // The actual gradient-descent update:
            Eigen::VectorXd trial = x - alpha * g;

            // Keep boundary positions exactly at their starting values.
            for (int vertex : pinned_vertices) {
                trial.segment<3>(3*vertex) =
                    current_components.segment<3>(3*vertex);
            }

            double energy_after = total_energy(trial);

            // Require sufficient decrease in energy.
            if (trial.allFinite() &&
                std::isfinite(energy_after) &&
                energy_after <=
                    energy_before - 1e-4 * alpha * g.squaredNorm())
            {
                x = trial;
                accepted = true;
                break;
            }

            // Retry with a smaller step.
            alpha *= 0.5;
        }

        if (!accepted) {
            throw std::runtime_error(
                "Could not find a step that decreases energy."
            );
        }
    }

    // Return the latest positions, even if not yet converged.
    return x;
}


//try gradient descent
//have codex vib code the visualization
//use houdini, every step of gradient descent, visualize as a frame
//have to code gravity
//minimize potential with gradient decent
//nonlinear gasuadel




int main(int argc, char* argv[])
{
    if (argc > 1) {
        if (std::string(argv[1]) != "--houdini" || argc > 3) {
            std::cerr << "Usage: hello_eigen [--houdini [output_directory]]\n";
            return 1;
        }
        try {
            export_houdini_demo(argc == 3 ? argv[2] : "output/houdini");
        } catch (const std::exception& error) {
            std::cerr << error.what() << '\n';
            return 1;
        }
        return 0;
    }
    int N = 2;
    int indicies = 4;
    double mew = 3.0;
    double lambda = 4.0;
    Eigen::VectorXi mesh(3*N);
    Eigen::VectorXd reference_components(3*indicies);
    Eigen::VectorXd current_components(3*indicies);
    mesh << 0, 1, 3,
             1, 2, 3;

    reference_components << 0, 0, 0,
                            1, 0, 0,
                            1, 1, 0,
                            0, 1, 0;

    current_components <<   0, 0, 0,
                            1, 0, 0,
                            1, 1, 0,
                            0, 1, 0;

    double PE1 = PE(N, mesh, reference_components, current_components, mew, lambda);

    std::cout << "\n Here is PE: " << PE1 << std::endl;

    Matrix32 F;
    F << 1, 2,
         3, 4,
         5, 6;

    // Eigen::MatrixXd F = Eigen::MatrixXd::Random(3, 2);

    double phi = psi(F, lambda, mew);
    Matrix32 D = D_psi(F, lambda , mew);

    test_psi();
    test_DD_psi();
    test_gradient();

    std::cout << "\n here is derivative of psi \n" << D <<std::endl;
    std::cout << "\n here is  psi \n" << phi <<std::endl;
    return 0;
}
