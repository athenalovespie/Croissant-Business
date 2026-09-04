#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <vector>


using Matrix32 = Eigen::Matrix<double, 3, 2>;
using Matrix22 = Eigen::Matrix<double, 2, 2>;
using Matrix33 = Eigen::Matrix<double, 3, 3>;

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
    std::cout << "Here is the matrix F:\n" << F << std::endl;
    
    auto [sigmas, U, V] = get_svd(F);
    double sigma_1 = sigmas(0);
    double sigma_2 = sigmas(1);

    std::cout << "\n here is U \n" << U << std::endl;
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

//takes in 3 points and then finds the area between them, these points are on the cartesian plane
//have X_i points , need all of the 

//make a long vector full of points X
//


double area(int t, Eigen::VectorXi& mesh, Eigen::VectorXd& reference_components, Eigen::VectorXd& current_components){
    int V1 = mesh(3*t);
    int V2 = mesh(3*t+1);
    int V3 = mesh(3*t+2);

    double x11 = current_components(3*V1);
    double x12 = current_components(3*V1+1);
    double x13 = current_components(3*V1+2);

    double x21 = current_components(3*V2);
    double x22 = current_components(3*V2+1);
    double x23 = current_components(3*V2+2);

    double x31 = current_components(3*V3);
    double x32 = current_components(3*V3+1);
    double x33 = current_components(3*V3+2);

    return 1/2*std::abs(x11*(x22 - x12) + x21*(x32 - x12) +x31*(x12 -x22));

    // return 1/2*std::abs(V1(0)*(V2(1)-V1(1)) + V2(0)*(V3(1) - V1(1)) + V3(0)*(V1(1) - V2(1)));
}

// double get_components(int a, int b, int t, Eigen::VectorXi& mesh, Eigen::VectorXd& current_components){
//     int vertex = mesh(3*t + (a-1));
//     double xab = current_components(3*vertex + (b-1));
//     return xab;
// }

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
    DS << (x11 - x33), (x21 - x31),
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



void test_psi() {
    Matrix32 F;
    F << 1, 2,
         3, 4,
         5, 6;

    Matrix32 direction;
    direction << 1, -1,
                 2,  1,
                -1,  2;
    direction.normalize();

    const double mu = 2.0;
    const double lambda = 3.0;
    const double energy = psi(F, lambda, mu);
    const Matrix32 gradient = D_psi(F, lambda, mu);

    const double directional_derivative = (gradient.array() * direction.array()).sum();

    constexpr int number_of_points = 100;
    Eigen::VectorXd x_axis(number_of_points);
    Eigen::VectorXd y_axis(number_of_points);

    int i = 0;
    for (int k = 10000; k>0; k-=100){
        double epk = 1.0 / static_cast<double>(k);
        double finite_difference =    //depends on psi
        (psi(F + epk * direction, lambda, mu)
        - psi(F - epk * direction, lambda, mu)) / (2.0 * epk);
        double error = std::abs((finite_difference)-(directional_derivative));
        x_axis[i] = epk;
        y_axis[i] = error;
        i++;
    }

    std::ofstream output_file("plot_data.csv");
    if (!output_file) {
        std::cerr << "Could not open plot_data.csv for writing.\n";
        return;
    }

    output_file << "x_axis,y_axis\n";
    output_file << std::setprecision(17);
    for (Eigen::Index row = 0; row < x_axis.size(); ++row) {
        output_file << x_axis(row) << ',' << y_axis(row) << '\n';
    }

    std::cout << "\nWrote " << x_axis.size()
              << " plotting points to plot_data.csv\n";

}

void test_DD_psi()
{
    Matrix32 F;
    F << 1, 2,
         3, 4,
         5, 6;

    Matrix32 delta_F;
    delta_F << 1, -1,
               2,  1,
              -1,  2;

    // Make it a unit direction
    delta_F.normalize();

    const double mu = 2.0;
    const double lambda = 3.0;

    // Your analytical second derivative applied to delta_F
    Matrix32 analytical =
        DD_psi(F, delta_F, lambda, mu);

    constexpr int number_of_points = 100;

    Eigen::VectorXd x_axis(number_of_points);
    Eigen::VectorXd y_axis(number_of_points);

    int i = 0;

    for (int k = 10000; k > 0; k -= 100)
    {
        double eps = 1.0 / static_cast<double>(k);

        // Finite difference of the FIRST derivative
        Matrix32 D_plus =
            D_psi(F + eps * delta_F, lambda, mu);

        Matrix32 D_minus =
            D_psi(F - eps * delta_F, lambda, mu);

        Matrix32 finite_difference =
            (D_plus - D_minus) / (2.0 * eps);

        // Frobenius norm of the difference
        double error =
            (finite_difference - analytical).norm();

        x_axis[i] = eps;
        y_axis[i] = error;

        i++;
    }

    std::ofstream output_file("DD_plot_data.csv");

    if (!output_file)
    {
        std::cerr
            << "Could not open DD_plot_data.csv for writing.\n";
        return;
    }

    output_file << "x_axis,y_axis\n";
    output_file << std::setprecision(17);

    for (Eigen::Index row = 0; row < x_axis.size(); ++row)
    {
        output_file
            << x_axis(row)
            << ','
            << y_axis(row)
            << '\n';
    }

    std::cout
        << "\nAnalytical DD_psi : delta_F:\n"
        << analytical
        << "\n";

    std::cout
        << "\nWrote "
        << x_axis.size()
        << " plotting points to DD_plot_data.csv\n";
}

int main()
{
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

    // Matrix32 F;
    // F << 1, 2,
    //      3, 4,
    //      5, 6;

    // double mew = 1.0;
    // double lambda = 1.0;
    // // Eigen::MatrixXd F = Eigen::MatrixXd::Random(3, 2);

    // double phi = psi(F, lambda, mew);
    // Matrix32 D = D_psi(F, lambda , mew);

    // test_psi();
    // test_DD_psi();

    // std::cout << "\n here is derivative of psi \n" << D <<std::endl;
    // std::cout << "\n here is  psi \n" << phi <<std::endl;
    // return 0;
}







// void test_psi() {
//     Matrix32 F;
//     F << 1, 2,
//          3, 4,
//          5, 6;

//     Matrix32 direction;
//     direction << 1, -1,
//                  2,  1,
//                 -1,  2;
//     direction.normalize();

//     const double mu = 2.0;
//     const double lambda = 3.0;
//     const double epsilon = 1e-6;

//     const double energy = psi(F, lambda, mu);
//     const Matrix32 gradient = D_psi(F, lambda, mu);
//     const double directional_derivative = (gradient.array() * direction.array()).sum(); //depends on d_psi
//     const double finite_difference =    //depends on psi
//         (psi(F + epsilon * direction, lambda, mu)
//        - psi(F - epsilon * direction, lambda, mu)) / (2.0 * epsilon);

//     std::cout << "psi(F): " << energy << "\n\n";
//     std::cout << "d_psi(F):\n" << gradient << "\n\n";
//     std::cout << "Gradient directional derivative: " << directional_derivative << "\n";
//     std::cout << "Finite-difference derivative: " << finite_difference << "\n";
//     std::cout << "Absolute error: "
//               << std::abs(directional_derivative - finite_difference) << "\n";
// }
