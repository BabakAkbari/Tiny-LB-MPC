/*
 * This file is part of the FreeRTOS port to Teensy boards.
 * Copyright (c) 2020 Timo Sandmann
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file    main.cpp
 * @brief   FreeRTOS example for Teensy boards
 * @author  Timo Sandmann
 * @date    17.05.2020
 */

#include "arduino_freertos.h"
#if _GCC_VERSION < 60100
#error "Compiler too old for std::thread support with FreeRTOS."
#endif

#include <thread>
#include <future>
#include <chrono>
#include <string>
#include <SD.h>
#include <ctime>
#include "tinympc/admm.hpp"
#include "tinympc/cf_teensy_interface.h"
#include "tinympc/gp.h"
#include "filter.h"

using namespace std::chrono_literals;

Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
Eigen::IOFormat SaveData(4, 0, ", ", "\n");
IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " << ", ";");




static lpf2pData cmdxfilt;
static lpf2pData cmdyfilt;
static lpf2pData cmdzfilt;

static Vector<float, 4> u_cmd;
static Eigen::VectorXd a_cmd(3);
static Eigen::VectorXd ac_a(3);
static Eigen::Vector3d ref_a;
static Eigen::Vector3d g;
static Eigen::VectorXd states(10);
static Eigen::VectorXd drags(2);
static Eigen::Vector3d omega;
static Eigen::Vector3d target_accel;
static Eigen::Vector3d target_jerk;
std::vector<double> rpy {0, 0, 0};
Eigen::VectorXd _c(3);
Eigen::VectorXd integral_pos_e(3);
Eigen::VectorXd integral_rpy_e(3);

Eigen::VectorXd P_COEFF_FOR(3);
Eigen::VectorXd I_COEFF_FOR(3);
Eigen::VectorXd D_COEFF_FOR(3);

Eigen::VectorXd P_COEFF_TOR(3); 
Eigen::VectorXd I_COEFF_TOR(3); 
Eigen::VectorXd D_COEFF_TOR(3); 
Eigen::VectorXd last_rpy_e(3);
Eigen::VectorXd last_pos_e(3);

int n = 20;
int d = 10;
int iter_conf = 3;

double sigma_eta_x = 0.001230735319648732;
double sigma_eta_y = 0.001230735319648732;
double sigma_eta_z = 0.00027692586973849045;
double sigma_n = 1;

// MatrixXd X(n, d);
VectorXd Yx(n);
VectorXd Yy(n);
VectorXd Yz(n);
MatrixXd M_inv_x = MatrixXd::Identity(d, d);
MatrixXd K_X_X_inverse_x;

MatrixXd M_inv_y = MatrixXd::Identity(d, d);
MatrixXd K_X_X_inverse_y;

MatrixXd M_inv_z = MatrixXd::Identity(d, d);
MatrixXd K_X_X_inverse_z;

Eigen::VectorXd diagonalValues_x(10);
Eigen::VectorXd diagonalValues_y(10);
Eigen::VectorXd diagonalValues_z(10);

std::vector<VectorXd> z_stars;
std::vector<VectorXd> m_z_x(10);
std::vector<VectorXd> m_z_y(10);
std::vector<VectorXd> m_z_z(10);
std::vector<MatrixXd> v_hat_x(10);
std::vector<MatrixXd> v_hat_y(10);
std::vector<MatrixXd> v_hat_z(10);
std::vector<MatrixXd> L_hat_x(10);
std::vector<MatrixXd> L_hat_y(10);
std::vector<MatrixXd> L_hat_z(10);

MatrixXd x_prev;
Eigen::MatrixXd Q_glob = MatrixXd::Zero(10, 10);
MatrixXd H_xz = MatrixXd::Zero(2, 10);
MatrixXd m_xz = MatrixXd::Zero(2, 11);
MatrixXd H_x = MatrixXd::Zero(1, 10);
MatrixXd m_x = MatrixXd::Zero(1, 11);

MatrixXd H_xyz = MatrixXd::Zero(3, 10);
MatrixXd m_xyz = MatrixXd::Zero(3, 11);
MatrixXd H_xy = MatrixXd::Zero(2, 10);
MatrixXd m_xy = MatrixXd::Zero(2, 11);

MatrixXd L_x = MatrixXd::Zero(11, 11);
MatrixXd L_x_r = MatrixXd::Zero(11, 10);
MatrixXd L_x_c = MatrixXd::Zero(11, 11);

MatrixXd L_y = MatrixXd::Zero(11, 11);
MatrixXd L_y_r = MatrixXd::Zero(11, 10);
MatrixXd L_y_c = MatrixXd::Zero(11, 11);

MatrixXd L_z = MatrixXd::Zero(11, 11);
MatrixXd L_z_r = MatrixXd::Zero(11, 10);
MatrixXd L_z_c = MatrixXd::Zero(11, 11);


MatrixXd ref = MatrixXd::Zero(10, 10);

Eigen::VectorXd z_bar(11);

Eigen::VectorXd _p(3);
Eigen::VectorXd _v(3);
Eigen::VectorXd _a(3);
Eigen::VectorXd _j(3);
Eigen::VectorXd _w(3);
MatrixXd _R = Eigen::MatrixXd::Identity(3, 3);
Eigen::VectorXd c(3);
Eigen::VectorXd tau(3);
Eigen::VectorXd ref_p(3);
Eigen::VectorXd ref_v(3);
Eigen::MatrixXd ref_R = Eigen::MatrixXd::Identity(3, 3);
double _dt_dyn = 0.002;

TinySolver *solver;


std::vector<MatrixXd> K(10);
std::vector<MatrixXd> P(10);
std::vector<MatrixXd> Quu_inv(10);
std::vector<MatrixXd> AmBKt(10);

static Eigen::VectorXd drag_est(2);
static Eigen::VectorXd drag_uncert(2);
static Eigen::VectorXd pos(6);
MatrixXd P_N = MatrixXd::Zero(10, 10);
static std::vector<Eigen::VectorXd> poses;
int flag = 0;
double DT = 0.01;
struct request req;
struct response res;
uint8_t req_buffer[60];
uint8_t res_buffer[60];
int last_connected = 0;

bool response_valid(struct response* res) {
    float* data = (float*)res;
    for(int i = 0; i < 15; i ++) {
        if ((data[i] != 0.0 && !data[i]) || data[i] > 100 || data[i] < -100) return false;
    }
    return true;
}

bool initializeSDCard() {
    if (!SD.begin(BUILTIN_SDCARD)) {
        Serial.println("SD card initialization failed!");
        return false;  // Initialization failed
    } else {
        Serial.println("SD card initialization done.");
        return true;   // Initialization successful
    }
}

void appendToCSV(const std::vector<Eigen::VectorXd> &data, const char *filename)
{
    File file = SD.open(filename, FILE_WRITE);

    if (!file)
    {
        Serial.println("Failed to open file for appending");
        return;
    }

    for (const auto &row : data)
    {
        for (int i = 0; i < row.size(); ++i)
        {
            file.print(row[i]);

            if (i < row.size() - 1)
            {
                file.print(",");
            }
        }
        file.println(); // Line break after each row
    }
    // file.println("\n");

    file.close();
}


void appendvecToCSV(const Eigen::VectorXd &data, const char* filename) {
    File file = SD.open(filename, FILE_WRITE);
    // Serial.println(filename);
    if (!file) {
        Serial.println("Error opening file");
        return;
    }

    for (int i = 0; i < data.size(); i++) {
        file.print(data[i], 6);  // You can specify the precision as needed
        if (i < data.size() - 1) {
            file.print(",");
        }
    }
    file.println();
    file.close();
}

void  deleteFileIfExists(const char* filename) {
    if (SD.exists(filename)) {
        if (SD.remove(filename)) {
            Serial.println("File deleted successfully");
        } else {
            Serial.println("Failed to delete the file");
        }
    } else {
        Serial.println("File does not exist");
    }
}

void readCSVFiles(const char* filePathX, const char* filePathY, 
                  Eigen::MatrixXd &X, Eigen::VectorXd &Yx, Eigen::VectorXd &Yz) {
    File fileX = SD.open(filePathX);
    File fileY = SD.open(filePathY);
    
    if (!fileX || !fileY) {
        Serial.println("Error opening files");
        return;
    }

    // Temporary storage to determine size and store data
    std::vector<std::vector<double>> X_data;
    std::vector<double> Yx_data, Yz_data;

    // Read and store data from X file
    while (fileX.available()) {
        String line = fileX.readStringUntil('\n');
        std::vector<double> row;
        int lastComma = -1;

        for (int i = 0; i < line.length(); i++) {
            if (line[i] == ',') {
                row.push_back(line.substring(lastComma + 1, i).toFloat());
                lastComma = i;
            }
        }
        row.push_back(line.substring(lastComma + 1).toFloat());
        X_data.push_back(row);
    }
    fileX.close();

    // Read and store data from Y file
    while (fileY.available()) {
        String line = fileY.readStringUntil('\n');
        int commaIdx = line.indexOf(',');

        Yx_data.push_back(line.substring(0, commaIdx).toFloat());
        Yz_data.push_back(line.substring(commaIdx + 1).toFloat());
    }
    fileY.close();

    // Convert vectors to Eigen structures
    int n = X_data.size();
    int d = X_data[0].size();

    X.resize(n, d);
    Yx.resize(n);
    Yz.resize(n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < d; ++j) {
            X(i, j) = X_data[i][j];
        }
        Yx(i) = Yx_data[i];
        Yz(i) = Yz_data[i];
    }
}

void readThetaFiles(const char* theta1Filename, const char* theta2Filename, 
                    double &sigma_eta_x, Eigen::VectorXd &diagonalValues_x,
                    double &sigma_eta_z, Eigen::VectorXd &diagonalValues_z) {
    // Open theta1.txt
    File theta1File = SD.open(theta1Filename);    if (theta1File) {
        String line = theta1File.readStringUntil('\n');
        std::vector<double> theta1ValuesVec;
        
        int lastComma = -1;
        for (int i = 0; i < line.length(); i++) {
            if (line[i] == ',') {
                theta1ValuesVec.push_back(line.substring(lastComma + 1, i).toFloat());
                lastComma = i;
            }
        }
        theta1ValuesVec.push_back(line.substring(lastComma + 1).toFloat());
        
        // Create Eigen::VectorXd from the vector
        Eigen::VectorXd theta1Values = Eigen::VectorXd::Map(theta1ValuesVec.data(), theta1ValuesVec.size());

        // Extract values
        sigma_eta_x = std::sqrt(theta1Values(0));
        diagonalValues_x = theta1Values.segment(1, 10);  // Assuming diagonalValues_x has 10 elements

        theta1File.close();
    } else {
        Serial.println("Unable to open theta1.txt");
        return;
    }

    // Open theta2.txt
    File theta2File = SD.open(theta2Filename);
    if (theta2File) {
        String line = theta2File.readStringUntil('\n');
        std::vector<double> theta2ValuesVec;

        int lastComma = -1;
        for (int i = 0; i < line.length(); i++) {
            if (line[i] == ',') {
                theta2ValuesVec.push_back(line.substring(lastComma + 1, i).toFloat());
                lastComma = i;
            }
        }
        theta2ValuesVec.push_back(line.substring(lastComma + 1).toFloat());

        // Create Eigen::VectorXd from the vector
        Eigen::VectorXd theta2Values = Eigen::VectorXd::Map(theta2ValuesVec.data(), theta2ValuesVec.size());

        // Extract values
        sigma_eta_z = std::sqrt(theta2Values(0));
        diagonalValues_z = theta2Values.segment(1, 10);  // Assuming diagonalValues_z has 10 elements

        theta2File.close();
    } else {
        Serial.println("Unable to open theta2.txt");
        return;
    }
}



void quadmodel(Eigen::VectorXd& _p, Eigen::VectorXd& _v, Eigen::VectorXd& _a, Eigen::VectorXd& _j, Eigen::VectorXd& _w, MatrixXd& _R, Eigen::VectorXd c,Eigen::VectorXd tau, double dt) {
    // Define Constatnts
    Eigen::VectorXd g(3);
    g << 0, 0, 9.8;
    double m = 1;
    Matrix3d J;
    J << 2.3951e-5, 0, 0,
         0, 2.3951e-5, 0,
         0, 0, 3.2347e-5;
    double env_density = 100;
    Eigen::VectorXd cd(3);
    cd << 0.5, 0.5, 0.5;
    Eigen::VectorXd cd_ref_area(3);
    cd_ref_area << 0.159, 0.159, 0.159;
    Eigen::VectorXd sq_drag(3);
    
    // Compute drag
    double sign_v_x = ((_R.transpose() * _v)[0] >= 0) ? 1.0 : -1.0;
    double sign_v_y = ((_R.transpose() * _v)[1] >= 0) ? 1.0 : -1.0;
    double sign_v_z = ((_R.transpose() * _v)[2] >= 0) ? 1.0 : -1.0;
    Eigen::VectorXd sign_v(3);
    sign_v << sign_v_x, sign_v_y, sign_v_z;
    sq_drag = 0.5 * env_density * sign_v.cwiseProduct(cd).cwiseProduct(cd_ref_area).cwiseProduct(_R.transpose() * _v).cwiseProduct(_R.transpose() * _v);
    Eigen::VectorXd drag = _R * sq_drag;

    MatrixXd D(3, 3);
    D << 5, 0, 0,
         0, 5, 0,
         0, 0, 5;
    drag = _R * D * _R.transpose() * _v;

    Eigen::VectorXd a = -g + _R * c - drag;
    // Dynamics Equations
    // Eigen::VectorXd a = (-m * g + _R * c - drag) / m;
    // Eigen::VectorXd a = -g +  _R * c;
    Eigen::VectorXd v = _v + a * dt;
    Eigen::VectorXd p = _p + v * dt;
    Vector3d _w3 = _w;
    Eigen::VectorXd w_dot = J.inverse() * (tau - _w3.cross(J * _w3)); 
    Eigen::VectorXd w = _w + w_dot * dt;
    MatrixXd skew_w(3, 3);
    skew_w <<  0, -w(2), w(1),
               w(2), 0, -w(0),
               -w(1), w(0), 0;
    MatrixXd R_dot(3, 3);
    R_dot = _R * skew_w;
    Eigen::Matrix3d rotationUpdate = (skew_w * dt).exp();
    MatrixXd R(3, 3);
    R = _R * rotationUpdate;
    Eigen::VectorXd ez(3);
    ez << 0, 0, 1;
    Eigen::VectorXd c_dot = (c - _c) / dt;
    Eigen::VectorXd j = R * c_dot + c(2) * R * skew_w * ez;
    if(p[2] < 0){
        p[2] = 0;
        v << 0, 0, 0;
        w << 0, 0, 0;
    }

    // Compute thrust with the next states
    sign_v_x = ((R.transpose() * v)[0] >= 0) ? 1.0 : -1.0;
    sign_v_y = ((R.transpose() * v)[1] >= 0) ? 1.0 : -1.0;
    sign_v_z = ((R.transpose() * v)[2] >= 0) ? 1.0 : -1.0;
    sign_v << sign_v_x, sign_v_y, sign_v_z;
    sq_drag = 0.5 * env_density * sign_v.cwiseProduct(cd).cwiseProduct(cd_ref_area).cwiseProduct(R.transpose() * v).cwiseProduct(R.transpose() * v);
    drag = R * sq_drag;
    drag = R * D * R.transpose() * v;
    drags << drag(0), drag(2); 
    
    // Update the states
    _p = p;
    _v = v;
    _a = a;
    _j = j;
    _w = w;
    _R = R;
    _c = c;
}

void computeomega(Eigen::Vector3d *omega, Eigen::Vector3d &target_x_ax, Eigen::Vector3d &target_y_ax, Eigen::Vector3d &target_accel, Eigen::Vector3d &target_jerk){
  double c = target_accel.norm();
  double term1 = -target_y_ax.transpose() * target_jerk;
  (*omega)(0) = (term1 / c) * (180.0 / M_PI);
  double term2 = target_x_ax.transpose() * target_jerk;
  (*omega)(1) = (term2 / c) * (180.0 / M_PI);
  (*omega)(2) = 0; 
}

Eigen::VectorXd rotationMatrixToEulerAngles(const Eigen::MatrixXd& R) {
    double roll, pitch, yaw;

    // Extract angles using atan2
    roll = std::atan2(R(2, 1), R(2, 2));
    pitch = std::atan2(-R(2, 0), std::sqrt(R(2, 1) * R(2, 1) + R(2, 2) * R(2, 2)));
    yaw = std::atan2(R(1, 0), R(0, 0));

    VectorXd rpy(3);
    rpy << roll, pitch, yaw;
    return rpy;
}

Eigen::Matrix3d rpyToRotationMatrix(const Eigen::Vector3d& rpy) {
    double roll = rpy[0];
    double pitch = rpy[1];
    double yaw = rpy[2];

    Eigen::AngleAxisd roll_rotation(roll, Eigen::Vector3d::UnitX());
    Eigen::AngleAxisd pitch_rotation(pitch, Eigen::Vector3d::UnitY());
    Eigen::AngleAxisd yaw_rotation(yaw, Eigen::Vector3d::UnitZ());

    Eigen::Quaterniond quaternion = yaw_rotation * pitch_rotation * roll_rotation;

    return quaternion.toRotationMatrix();
}

double squaredExponentialKernel(const VectorXd &zi, const VectorXd &zj, const MatrixXd &M_inv,
                                double sigma_eta, double sigma_n, int i, int j)
{
    VectorXd diff = zi - zj;
    // double distanceSquared = diff.dot(M_inv * diff);
    double distanceSquared = diff.transpose() * M_inv * M_inv * diff;
    // std::cout << "diff " << distanceSquared;
    double kernelValue = sigma_eta * sigma_eta * exp(-0.5 * distanceSquared);

    if (i == j)
    {
        kernelValue += sigma_n * sigma_n;
    }

    return kernelValue;
}

VectorXd computem(const VectorXd &z_star, const MatrixXd &X, const VectorXd &Y, const MatrixXd &M_inv, MatrixXd K_X_X_inverse, double sigma_eta, double sigma_n)
{

    int N = X.rows();
    int n = z_star.size();

    VectorXd mz(n + 1);

    // first matrix in term 2
    MatrixXd term1_1 = MatrixXd::Zero(n + 1, N);
    VectorXd k_zstar_X(N);

    for (int i = 0; i < N; i++)
    {
        k_zstar_X(i) = squaredExponentialKernel(z_star, X.row(i), M_inv, sigma_eta, sigma_n, -1, i);
    }

    term1_1.block(0, 0, 1, N) = k_zstar_X.transpose();

    MatrixXd k_zstar_X10 = MatrixXd::Zero(n, N);

    for (int i = 0; i < N; i++)
    {
        // k_zstar_X10.col(i) = (M_inv * M_inv) * (X.row(i).transpose() - z_star) * squaredExponentialKernel(z_star, X.row(i), M_inv, sigma_eta, sigma_n, -1, i);
        k_zstar_X10.col(i) = (M_inv * M_inv) * (X.row(i).transpose() - z_star) * k_zstar_X(i);
    }

    term1_1.block(1, 0, n, N) = k_zstar_X10;

    // seccond matrix in term 2

    // MatrixXd K_X_X = MatrixXd::Zero(N, N); // K(X, X)
    // for (int i = 0; i < N; ++i)
    // {
    //     for (int j = 0; j < N; ++j)
    //     {
    //         K_X_X(i, j) = squaredExponentialKernel(X.row(i), X.row(j), M_inv, sigma_eta, 0.0, i, j);
    //     }
    // }

    // MatrixXd K_X_X_inverse = (K_X_X + sigma_eta * sigma_eta * MatrixXd::Identity(N, N)).inverse();

    mz = term1_1 * K_X_X_inverse * Y;

    return mz;
}

MatrixXd computeVhat(const MatrixXd &z_star, const MatrixXd &X, const MatrixXd &M_inv, MatrixXd K_X_X_inverse,
                     double sigma_eta, double sigma_n) 
{
    int N = X.rows();
    int n = z_star.size();

    MatrixXd v_hat = MatrixXd::Zero(n + 1, n + 1);

    // term 1

    MatrixXd term1 = MatrixXd::Zero(n + 1, n + 1);
    // Compute k(z*,z*) and set it in the matrix term1
    term1(0, 0) = sigma_eta * sigma_eta + sigma_n * sigma_n;

    // Compute K^{(0,1)}(z*,z*) and set it in the matrix term1
    VectorXd K_01(n);
    K_01.setZero();
    // for (int i = 0; i < n; ++i)
    // {
    //     K_01(i) = 0;
    // }
    term1.block(0, 1, 1, n) = K_01.transpose();

    // Compute K^{(1,0)}(z*,z*) and set it in the matrix term1
    VectorXd K_10(n);
    K_10.setZero();
    // for (int i = 0; i < n; ++i)
    // {
    //     K_10(i) = 0;
    // }
    term1.block(1, 0, n, 1) = K_10;

    // Compute K^{(1,1)}(z*,z*) and set it in the matrix term1
    MatrixXd K_11 = (sigma_eta * sigma_eta) * (M_inv * M_inv);
    term1.block(1, 1, n, n) = K_11;

    // end of term1

    // term 2
    MatrixXd term2 = MatrixXd::Zero(n + 1, n + 1);

    // first matrix in term 2
    MatrixXd term2_1 = MatrixXd::Zero(n + 1, N);
    VectorXd k_zstar_X(N);

    for (int i = 0; i < N; i++)
    {
        k_zstar_X(i) = squaredExponentialKernel(z_star, X.row(i), M_inv, sigma_eta, sigma_n, -1, i);
    }
    term2_1.block(0, 0, 1, N) = k_zstar_X.transpose();

    MatrixXd k_zstar_X10 = MatrixXd::Zero(n, N);

    for (int i = 0; i < N; i++)
    {
        // k_zstar_X10.col(i) = (M_inv * M_inv) * (X.row(i).transpose() - z_star) * squaredExponentialKernel(z_star, X.row(i), M_inv, sigma_eta, sigma_n, -1, i);
        k_zstar_X10.col(i) = (M_inv * M_inv) * (X.row(i).transpose() - z_star) * k_zstar_X(i);
    }

    term2_1.block(1, 0, n, N) = k_zstar_X10;

    // seccond matrix in term 2

    // MatrixXd K_X_X = MatrixXd::Zero(N, N); // K(X, X)
    // for (int i = 0; i < N; ++i)
    // {
    //     for (int j = 0; j < N; ++j)
    //     {
    //         K_X_X(i, j) = squaredExponentialKernel(X.row(i), X.row(j), M_inv, sigma_eta, 0.0, i, j);
    //     }
    // }

    // MatrixXd K_X_X_inverse = (K_X_X + sigma_eta * sigma_eta * MatrixXd::Identity(N, N)).inverse();

    // third matrix in term 2
    MatrixXd term2_3 = MatrixXd::Zero(N, n + 1);
    VectorXd k_X_zstar(N);

    // for (int i = 0; i < N; i++)
    // {
    //     k_X_zstar(i) = squaredExponentialKernel(X.row(i), z_star, M_inv, sigma_eta, sigma_n, -1, i);
    //     // std::cout<< "k_xzstar " << k_X_zstar(i) << std::endl;
    // }
    // term2_3.block(0, 0, N, 1) = k_X_zstar;
    term2_3.block(0, 0, N, 1) = k_zstar_X;
    // MatrixXd k_X_zstar01 = MatrixXd::Zero(N, n);

    // for (int i = 0; i < N; i++)
    // {
    //     // k_X_zstar01.row(i) = (M_inv * M_inv) * (X.row(i).transpose() - z_star) * squaredExponentialKernel(X.row(i), z_star, M_inv, sigma_eta, sigma_n, -1, i);
    //     k_X_zstar01.row(i) = (M_inv * M_inv) * (X.row(i).transpose() - z_star) * k_zstar_X(i);
    // }

    // term2_3.block(0, 1, N, n) = k_X_zstar01;

    term2_3.block(0, 1, N, n) = k_zstar_X10.transpose();
    term2 = term2_1 * K_X_X_inverse * term2_3;

    // end of term2

    v_hat = term1 - term2;
    return v_hat;
}

void PIDPositionControl(Eigen::VectorXd curr_p, Eigen::VectorXd curr_v, Eigen::VectorXd curr_a, Eigen::MatrixXd curr_R, Eigen::VectorXd ref_p, Eigen::VectorXd ref_v, Eigen::VectorXd ref_a, Eigen::MatrixXd ref_R, Eigen::VectorXd& target_thrust, MatrixXd& target_rotation, double dt) {

    Eigen::VectorXd g(3);
    g << 0, 0, 9.8;
    Eigen::VectorXd pos_e = ref_p - curr_p;
    Eigen::VectorXd vel_e =  ref_v - curr_v;
    Eigen::VectorXd d_pos_e = (pos_e - last_pos_e) / dt;
    last_pos_e = pos_e;
    integral_pos_e = integral_pos_e + pos_e * dt;
    integral_pos_e[0] = std::min(std::max(integral_pos_e[0], -2.), 2.);
    integral_pos_e[1] = std::min(std::max(integral_pos_e[1], -2.), 2.);
    integral_pos_e[2] = std::min(std::max(integral_pos_e[2], -2.), 2.);

    Eigen::VectorXd target_force = g
                                + P_COEFF_FOR.cwiseProduct(pos_e) 
                                + I_COEFF_FOR.cwiseProduct(integral_pos_e) 
                                + D_COEFF_FOR.cwiseProduct(vel_e);
    // target_force = ref_a + g;

    Eigen::VectorXd target_rpy(3);
    Eigen::Vector3d target_z_ax = target_force.normalized();
    Eigen::Vector3d target_x_c(std::cos(0), std::sin(0), 0);
    Eigen::Vector3d target_y_ax = target_z_ax.cross(target_x_c).normalized();
    Eigen::Vector3d target_x_ax = target_y_ax.cross(target_z_ax);
    target_rotation << target_x_ax, target_y_ax, target_z_ax;
    
    target_thrust = curr_R.transpose() * target_force;
    // target_thrust[2] = std::min(std::max(target_thrust[2], 0.0), 30.0);
    target_thrust << 0, 0, target_thrust[2];
}

void PIDAttitudeControl(Eigen::MatrixXd curr_R, Eigen::MatrixXd ref_R, Eigen::VectorXd& target_torques, double dt) {

    Eigen::VectorXd cur_rpy = rotationMatrixToEulerAngles(curr_R);
    Eigen::VectorXd target_rpy = rotationMatrixToEulerAngles(ref_R);
    Eigen::VectorXd rpy_e = target_rpy - cur_rpy;
    if (rpy_e(2) > M_PI)
        rpy_e(2) = rpy_e(2) - 2 * M_PI;
     if (rpy_e(2) < -M_PI)
        rpy_e(2) = rpy_e(2) + 2 * M_PI;
    Eigen::VectorXd d_rpy_e = (rpy_e - last_rpy_e) / dt;
    last_rpy_e = rpy_e;
    integral_rpy_e = integral_rpy_e + rpy_e * dt;

    target_torques = P_COEFF_TOR.cwiseProduct(rpy_e) 
                    + I_COEFF_TOR.cwiseProduct(integral_rpy_e)
                    + D_COEFF_TOR.cwiseProduct(d_rpy_e);
}

void PIDrateControl(Eigen::VectorXd curr_w, Eigen::VectorXd ref_w, Eigen::VectorXd& target_torques, double dt) {

    Eigen::VectorXd cur_rpy = curr_w;
    Eigen::VectorXd target_rpy = ref_w;
    Eigen::VectorXd rpy_e = target_rpy - cur_rpy;

    Eigen::VectorXd d_rpy_e = (rpy_e - last_rpy_e) / dt;
    last_rpy_e = rpy_e;
    integral_rpy_e = integral_rpy_e + rpy_e * dt;

    target_torques = P_COEFF_TOR.cwiseProduct(rpy_e) 
                    + I_COEFF_TOR.cwiseProduct(integral_rpy_e)
                    + D_COEFF_TOR.cwiseProduct(d_rpy_e);
    Matrix3d J;
    J << 2.3951e-5, 0, 0,
         0, 2.3951e-5, 0,
         0, 0, 3.2347e-5;
    target_torques = J * target_torques;
}

void dynamicsLoop() {
    TickType_t xLastWakeTime = xTaskGetTickCount();
    const TickType_t xFrequency = 2 / portTICK_PERIOD_MS;
    // 500 Hz Dynamics Loop
    while (true) {
        unsigned long startTime = millis();
        // PIDAttitudeControl(_R, target_rotation, tau, _dt_dyn);
        states << _p(0), _p(1), _p(2), _v(0), _v(1), _v(2), _a(0), _a(1), _a(2), 0;
        PIDrateControl(_w, omega, tau, _dt_dyn); 
        quadmodel(_p, _v, _a, _j, _w, _R, c, tau, _dt_dyn);
        vTaskDelayUntil(&xLastWakeTime, xFrequency);
        unsigned long endTime = millis();
        unsigned long duration = endTime - startTime;
    }
}

void gpLoop() {
    // if (solver->settings->en_state_bound == 0)
    //     ::vTaskSuspend(nullptr);

    TickType_t xLastWakeTime = xTaskGetTickCount();
    const TickType_t xFrequency = 100 / portTICK_PERIOD_MS;

    MatrixXd Q1;
    MatrixXd R1 = solver->work->R.array().matrix().asDiagonal();
    while (true) {
        unsigned long startTime = millis(); 
        
        // for (int j = 0; j < solver->work->N; j++) {
        //     if (j == solver->work->N - 1)
        //         solver->work->z_star.col(j) = solver->work->x.col(j);
        //     else
        //         solver->work->z_star.col(j) = solver->work->x.col(j+1);
        // }

        for (int j = 0; j < solver->work->N; j++) {
            solver->work->z_star.col(j) = ref.col(j);
        }

        int j = 9;
        v_hat_x[j] = computeVhat(solver->work->z_star.col(j), X_x, M_inv_x, K_X_X_inverse_x, sigma_eta_x, sigma_n);
        LLT<MatrixXd> llt_x(v_hat_x[j]);
        L_hat_x[j] = llt_x.matrixL();

        v_hat_y[j] = computeVhat(solver->work->z_star.col(j), X_x, M_inv_y, K_X_X_inverse_y, sigma_eta_y, sigma_n);
        LLT<MatrixXd> llt_y(v_hat_y[j]);
        L_hat_y[j] = llt_y.matrixL();
        
        v_hat_z[j] = computeVhat(solver->work->z_star.col(j), X_x, M_inv_z, K_X_X_inverse_z, sigma_eta_z, sigma_n);
        LLT<MatrixXd> llt_z(v_hat_z[j]);
        L_hat_z[j] = llt_z.matrixL();
        
        m_z_x[j] = computem(solver->work->z_star.col(j), X_x, Yx, M_inv_x, K_X_X_inverse_x, sigma_eta_x, sigma_n);
        m_z_y[j] = computem(solver->work->z_star.col(j), X_x, Yy, M_inv_y, K_X_X_inverse_y, sigma_eta_y, sigma_n);
        m_z_z[j] = computem(solver->work->z_star.col(j), X_x, Yz, M_inv_z, K_X_X_inverse_z, sigma_eta_z, sigma_n);
        // m_z_x[j].setZero();
        // m_z_y[j].setZero();
        // // m_z_z[j].setZero();
        // L_hat_x[j].setZero();
        // L_hat_y[j].setZero();  
        // L_hat_z[j].setZero();   
        if (solver->settings->en_state_bound == 1) {
            H_xyz.row(0) = m_z_x[j].segment(1, 10);
            H_xyz(0, 6) = m_z_x[j](7) + 1;
            H_xyz.row(1) = m_z_y[j].segment(1, 10);
            H_xyz(1, 7) = m_z_y[j](8) + 1;
            H_xyz.row(2) = m_z_z[j].segment(1, 10);
            H_xyz(2, 8) = m_z_z[j](9) + 1;

            m_xyz.row(0) = m_z_x[j];
            m_xyz.row(1) = m_z_y[j];
            m_xyz.row(2) = m_z_z[j];
            m_xyz(2, 0) = m_z_z[j](0) + 9.81;

            z_bar(0) = 1;
            z_bar.tail(solver->work->nx) = -solver->work->z_star.col(j);
            // z_bar.tail(solver->work->nx).setZero();

            solver->work->socx[0].A = H_xyz;
            solver->work->socx[0].b = m_xyz * z_bar;

            // L_x.setZero();
            L_x_r = L_hat_x[j].block(0, 1, 11, 10);
            solver->work->socx[1].A = L_x_r;
            solver->work->socx[1].b = L_hat_x[j] * z_bar;

            // L_z.setZero();
            L_z_r = L_hat_z[j].block(0, 1, 11, 10);
            solver->work->socx[2].A = L_z_r;
            solver->work->socx[2].b = L_hat_z[j] * z_bar;

            L_y_r = L_hat_y[j].block(0, 1, 11, 10);
            solver->work->socx[4].A = L_y_r;
            solver->work->socx[4].b = L_hat_y[j] * z_bar;

            H_xy.row(0) = H_xyz.row(0);
            H_xy.row(1) = H_xyz.row(1);
            m_xy.row(0) = m_xyz.row(0);
            m_xy.row(1) = m_xyz.row(1);

            solver->work->socx[3].A = H_xy;
            solver->work->socx[3].b = m_xy * z_bar;

            solver->work->H_z.col(j) = H_xyz.row(2);
            solver->work->m_hat_z.col(j) = m_xyz.row(2) * z_bar;

            solver->work->SC[j] = solver->work->socx;
            
            // SOC 0
            MatrixXd rhoI = solver->cache->rhox * MatrixXd::Identity(solver->work->socx[0].A.transpose().cols(), solver->work->socx[0].A.transpose().cols());
            MatrixXd ATrhOIA = solver->work->socx[0].A.transpose() * rhoI * solver->work->socx[0].A;
            ATrhOIA += solver->work->socx[0].A.transpose() * rhoI * solver->work->socx[0].A;
            ATrhOIA += solver->work->socx[0].A.transpose() * rhoI * solver->work->socx[0].A;
            
            // SOC 1
            rhoI = solver->cache->rhox * MatrixXd::Identity(solver->work->socx[1].A.transpose().cols(), solver->work->socx[1].A.transpose().cols());
            ATrhOIA += solver->work->socx[1].A.transpose() * rhoI * solver->work->socx[1].A;
            
            // SOC 2
            rhoI = solver->cache->rhox * MatrixXd::Identity(solver->work->socx[2].A.transpose().cols(), solver->work->socx[2].A.transpose().cols());
            ATrhOIA += solver->work->socx[2].A.transpose() * rhoI * solver->work->socx[2].A;
            
            // SOC 3
            rhoI = solver->cache->rhox * MatrixXd::Identity(solver->work->socx[3].A.transpose().cols(), solver->work->socx[3].A.transpose().cols());
            ATrhOIA += solver->work->socx[3].A.transpose() * rhoI * solver->work->socx[3].A;
            ATrhOIA += solver->work->socx[3].A.transpose() * rhoI * solver->work->socx[3].A;
            
            // SOC 4
            rhoI = solver->cache->rhox * MatrixXd::Identity(solver->work->socx[4].A.transpose().cols(), solver->work->socx[4].A.transpose().cols());
            ATrhOIA += solver->work->socx[4].A.transpose() * rhoI * solver->work->socx[4].A;
            
            // solver->work->Q[j]  << 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1000;
            // solver->work->Q[j]  << 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1000;

            solver->work->Q[j] = Q_glob;

            solver->work->Q[j] += ATrhOIA;
            Q1 = solver->work->Q[j];
            R1 = solver->work->R.array().matrix().asDiagonal();
            P[j] = P_N;
            solver->cache->Pinf[j] = P[j]; 
        }
        if (solver->settings->en_state_bound == 5) {
            // SOC 5
            MatrixXd rhoI = solver->cache->rhox * MatrixXd::Identity(solver->work->socx[5].A.transpose().cols(), solver->work->socx[5].A.transpose().cols());
            MatrixXd ATrhOIA = solver->work->socx[5].A.transpose() * rhoI * solver->work->socx[5].A;
            rhoI = solver->cache->rhox * MatrixXd::Identity(solver->work->socx[5].c.cols(), solver->work->socx[5].c.cols());
            ATrhOIA += solver->work->socx[5].c * rhoI * solver->work->socx[5].c.transpose();
            
            // solver->work->Q[j]  << 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1000;
            solver->work->Q[j] = Q_glob;

            solver->work->Q[j] += ATrhOIA;

            Q1 = solver->work->Q[j];
            R1 = solver->work->R.array().matrix().asDiagonal();
            P[j] = P_N;
            solver->cache->Pinf[j] = P[j]; 
        }       
        for (int j = solver->work->N - 2; j >= 0; j--) {
            // // std::cout << "z_star" <<solver->work->z_star.col(j)(0) << std::endl;
            // v_hat_x[j] = computeVhat(solver->work->z_star.col(j), X_x, M_inv_x, K_X_X_inverse_x, sigma_eta_x, sigma_n);
            // LLT<MatrixXd> llt_x(v_hat_x[j]);
            // L_hat_x[j] = llt_x.matrixL();

            // v_hat_y[j] = computeVhat(solver->work->z_star.col(j), X_x, M_inv_y, K_X_X_inverse_y, sigma_eta_y, sigma_n);
            // LLT<MatrixXd> llt_y(v_hat_y[j]);
            // L_hat_y[j] = llt_y.matrixL();
            
            // v_hat_z[j] = computeVhat(solver->work->z_star.col(j), X_x, M_inv_z, K_X_X_inverse_z, sigma_eta_z, sigma_n);
            // LLT<MatrixXd> llt_z(v_hat_z[j]);
            // L_hat_z[j] = llt_z.matrixL();
            
            // m_z_x[j] = computem(solver->work->z_star.col(j), X_x, Yx, M_inv_x, K_X_X_inverse_x, sigma_eta_x, sigma_n);
            // m_z_y[j] = computem(solver->work->z_star.col(j), X_x, Yy, M_inv_y, K_X_X_inverse_y, sigma_eta_y, sigma_n);
            // m_z_z[j] = computem(solver->work->z_star.col(j), X_x, Yz, M_inv_z, K_X_X_inverse_z, sigma_eta_z, sigma_n);

            L_hat_x[j] = L_hat_x[9];
            L_hat_y[j] = L_hat_y[9];
            L_hat_z[j] = L_hat_z[9];
            m_z_x[j] = m_z_x[9];
            m_z_y[j] = m_z_y[9];   
            m_z_z[j] = m_z_z[9]; 
 
        if (solver->settings->en_state_bound == 1) {
            H_xyz.row(0) = m_z_x[j].segment(1, 10);
            H_xyz(0, 6) = m_z_x[j](7) + 1;
            H_xyz.row(1) = m_z_y[j].segment(1, 10);
            H_xyz(1, 7) = m_z_y[j](8) + 1;
            H_xyz.row(2) = m_z_z[j].segment(1, 10);
            H_xyz(2, 8) = m_z_z[j](9) + 1;

            m_xyz.row(0) = m_z_x[j];
            m_xyz.row(1) = m_z_y[j];
            m_xyz.row(2) = m_z_z[j];
            m_xyz(2, 0) = m_z_z[j](0) + 9.81;

            z_bar(0) = 1;
            z_bar.tail(solver->work->nx) = -solver->work->z_star.col(j);
            // z_bar.tail(solver->work->nx).setZero();

            solver->work->socx[0].A = H_xyz;
            solver->work->socx[0].b = m_xyz * z_bar;

            // L_x.setZero();
            L_x_r = L_hat_x[j].block(0, 1, 11, 10);
            solver->work->socx[1].A = L_x_r;
            solver->work->socx[1].b = L_hat_x[j] * z_bar;

            // L_z.setZero();
            L_z_r = L_hat_z[j].block(0, 1, 11, 10);
            solver->work->socx[2].A = L_z_r;
            solver->work->socx[2].b = L_hat_z[j] * z_bar;

            L_y_r = L_hat_y[j].block(0, 1, 11, 10);
            solver->work->socx[4].A = L_y_r;
            solver->work->socx[4].b = L_hat_y[j] * z_bar;

            H_xy.row(0) = H_xyz.row(0);
            H_xy.row(1) = H_xyz.row(1);
            m_xy.row(0) = m_xyz.row(0);
            m_xy.row(1) = m_xyz.row(1);

            solver->work->socx[3].A = H_xy;
            solver->work->socx[3].b = m_xy * z_bar;

            solver->work->H_z.col(j) = H_xyz.row(2);
            solver->work->m_hat_z.col(j) = m_xyz.row(2) * z_bar;

            solver->work->SC[j] = solver->work->socx;
            

            // SOC 0
            MatrixXd rhoI = solver->cache->rhox * MatrixXd::Identity(solver->work->socx[0].A.transpose().cols(), solver->work->socx[0].A.transpose().cols());
            MatrixXd ATrhOIA = solver->work->socx[0].A.transpose() * rhoI * solver->work->socx[0].A;
            ATrhOIA += solver->work->socx[0].A.transpose() * rhoI * solver->work->socx[0].A;
            ATrhOIA += solver->work->socx[0].A.transpose() * rhoI * solver->work->socx[0].A;
            
            // SOC 1
            rhoI = solver->cache->rhox * MatrixXd::Identity(solver->work->socx[1].A.transpose().cols(), solver->work->socx[1].A.transpose().cols());
            ATrhOIA += solver->work->socx[1].A.transpose() * rhoI * solver->work->socx[1].A;
            
            // SOC 2
            rhoI = solver->cache->rhox * MatrixXd::Identity(solver->work->socx[2].A.transpose().cols(), solver->work->socx[2].A.transpose().cols());
            ATrhOIA += solver->work->socx[2].A.transpose() * rhoI * solver->work->socx[2].A;
            
            // SOC 3
            rhoI = solver->cache->rhox * MatrixXd::Identity(solver->work->socx[3].A.transpose().cols(), solver->work->socx[3].A.transpose().cols());
            ATrhOIA += solver->work->socx[3].A.transpose() * rhoI * solver->work->socx[3].A;
            ATrhOIA += solver->work->socx[3].A.transpose() * rhoI * solver->work->socx[3].A;

            // SOC 4
            rhoI = solver->cache->rhox * MatrixXd::Identity(solver->work->socx[4].A.transpose().cols(), solver->work->socx[4].A.transpose().cols());
            ATrhOIA += solver->work->socx[4].A.transpose() * rhoI * solver->work->socx[4].A;
            
            solver->work->Q[j] = Q_glob;

            solver->work->Q[j] += ATrhOIA;

            Q1 = solver->work->Q[j];

            K[j] = (R1 + solver->work->Bdyn.transpose() * P[j+1] * solver->work->Bdyn).inverse() * solver->work->Bdyn.transpose() * P[j+1] * solver->work->Adyn;
            P[j] = Q1 + solver->work->Adyn.transpose() * P[j+1] * (solver->work->Adyn - solver->work->Bdyn * K[j]);
            Quu_inv[j] = (R1 + solver->work->Bdyn.transpose() * P[j] * solver->work->Bdyn).inverse();
            AmBKt[j] = (solver->work->Adyn - solver->work->Bdyn * K[j]).transpose();
            
            // Compute cached matrices
            solver->cache->Kinf[j] = K[j];
            solver->cache->Pinf[j] = P[j];
            solver->cache->Quu_inv[j] = Quu_inv[j];
            solver->cache->AmBKt[j] = AmBKt[j];
            // // std::cout << "kinf is " << j << " : "<< solver->cache->Kinf[j] << std::endl; 
        }
        if(solver->settings->en_state_bound == 5) {
            // SOC 5
            MatrixXd rhoI = solver->cache->rhox * MatrixXd::Identity(solver->work->socx[5].A.transpose().cols(), solver->work->socx[5].A.transpose().cols());
            MatrixXd ATrhOIA = solver->work->socx[5].A.transpose() * rhoI * solver->work->socx[5].A;
            rhoI = solver->cache->rhox * MatrixXd::Identity(solver->work->socx[5].c.cols(), solver->work->socx[5].c.cols());
            ATrhOIA += solver->work->socx[5].c * rhoI * solver->work->socx[5].c.transpose();
            
            solver->work->Q[j] = Q_glob;

            solver->work->Q[j] += ATrhOIA;

            Q1 = solver->work->Q[j];

            K[j] = (R1 + solver->work->Bdyn.transpose() * P[j+1] * solver->work->Bdyn).inverse() * solver->work->Bdyn.transpose() * P[j+1] * solver->work->Adyn;
            P[j] = Q1 + solver->work->Adyn.transpose() * P[j+1] * (solver->work->Adyn - solver->work->Bdyn * K[j]);
            Quu_inv[j] = (R1 + solver->work->Bdyn.transpose() * P[j] * solver->work->Bdyn).inverse();
            AmBKt[j] = (solver->work->Adyn - solver->work->Bdyn * K[j]).transpose();
            
            // Compute cached matrices
            solver->cache->Kinf[j] = K[j];
            solver->cache->Pinf[j] = P[j];
            solver->cache->Quu_inv[j] = Quu_inv[j];
            solver->cache->AmBKt[j] = AmBKt[j];
        }
        }
        
        vTaskDelayUntil(&xLastWakeTime, pdMS_TO_TICKS(100));
        unsigned long endTime = millis();
        unsigned long duration = endTime - startTime;
    }
}

static void GPTask(void*) {
    gpLoop();
}

// static void DynamicsTask(void*) {
//     dynamicsLoop();
//     ::vTaskSuspend(nullptr);
// }

static void ControllerTask(void*) {
    // initializeSDCard(); 
    // deleteFileIfExists("pos.csv");
    // deleteFileIfExists("acc.csv");
    // deleteFileIfExists("ac_acc.csv");
    // deleteFileIfExists("states_train_gp.csv");
    // deleteFileIfExists("drags_train_gp.csv");
    // deleteFileIfExists("gptime.csv");
    // deleteFileIfExists("drag_vis.csv");
    // deleteFileIfExists("drag_est_vis.csv");
    // deleteFileIfExists("drag_uncert_vis.csv");
    // readThetaFiles("theta1.txt", "theta2.txt", sigma_eta_x, diagonalValues_x, sigma_eta_z, diagonalValues_z);

    sigma_eta_x = 0.1276430603884467;
    diagonalValues_x << 67368.76294399201,1.03599399490866,0.021337656615263847,10000000.000000006,1.5851640287877595,0.5468981234686067,1.1047747062668405,22086.34705927832,10000000.000000006,1.1644513400582364;
    sigma_eta_y = 0.017096395559383637;
    diagonalValues_y << 114.93186292240473,1.7352097795205164,0.050782774297603944,1431.427553355189,0.08463056198799858,86120.1694679899,214575.90236294555,20978.027668605093,0.21463707834634063,0.0016260381945378494;
    sigma_eta_z = 2.73483850811995;
    diagonalValues_z << 1.948445569299734,22.030778638193745,2670.542772022275,0.4550390196063133,28.23109551360207,0.6107564011088891,1.0477686453392376,10000.00000000001,2.000505035628158,135.60569995077;


    Yx = Y_x.col(0);
    Yy = Y_x.col(1);
    Yz = Y_x.col(2);



    
    // readCSVFiles("X_train_sampled.csv", "Y_train_sampled.csv", X, Yx, Yz);

    // Set the diagonal of M_inv to the specified values
    M_inv_x.diagonal() = diagonalValues_x;

    M_inv_x = M_inv_x.inverse();

    // std::cout << M_inv_x << std::endl;

    MatrixXd K_X_X_x = MatrixXd::Zero(n, n); // K(X, X)
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            K_X_X_x(i, j) = squaredExponentialKernel(X_x.row(i), X_x.row(j), M_inv_x, sigma_eta_x, 0.0, i, j);
        }
    }

    K_X_X_inverse_x = (K_X_X_x + sigma_n * sigma_n * MatrixXd::Identity(n, n)).inverse();

    M_inv_y.diagonal() = diagonalValues_y;

    M_inv_y = M_inv_y.inverse();

    // std::cout << M_inv_x << std::endl;

    MatrixXd K_X_X_y = MatrixXd::Zero(n, n); // K(X, X)
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            K_X_X_y(i, j) = squaredExponentialKernel(X_y.row(i), X_y.row(j), M_inv_y, sigma_eta_y, 0.0, i, j);
        }
    }

    K_X_X_inverse_y = (K_X_X_y + sigma_n * sigma_n * MatrixXd::Identity(n, n)).inverse();


    // Set the diagonal of M_inv to the specified values
    M_inv_z.diagonal() = diagonalValues_z;

    M_inv_z = M_inv_z.inverse();

    // std::cout << M_inv_z << std::endl;

    MatrixXd K_X_X_z = MatrixXd::Zero(n, n); // K(X, X)
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            K_X_X_z(i, j) = squaredExponentialKernel(X_y.row(i), X_y.row(j), M_inv_z, sigma_eta_z, 0.0, i, j);
        }
    }

    K_X_X_inverse_z = (K_X_X_z + sigma_n * sigma_n * MatrixXd::Identity(n, n)).inverse();
    
    int counter = 0;
    TickType_t xLastWakeTime = xTaskGetTickCount();
    const TickType_t xFrequency = 10 / portTICK_PERIOD_MS; // 10 ms
    Serial1.begin(2250000);
    while(!Serial1);
    Serial1.setTimeout(1000 * DT);

    // double freq = 3;
    // double period = ((5 * M_PI) / 2.0) / freq;

    // double freq = 2;
    // double period = ((5 * M_PI) / 2.0) / freq;
    // int end = (period / DT) + 1;
    
    // double freq = 1;
    // double period = (2 * M_PI) / freq;
    double freq = 1;
    double period = (6 * M_PI) / freq;
    int end = (period / DT) + 1;
    int sweepperiod = end / 6; // Position (sweep)
    double x_increment = 0.5 / 6;

    double theta = 45;  // 30 degrees in radians
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);
    
    // 100 Hz MPC Loop;
    while(true) {
        vTaskDelayUntil(&xLastWakeTime, pdMS_TO_TICKS(10));
        
        // Create buffer from request
        memcpy(req_buffer, &req, 60);
        // Sometimes there is junk data in the buffer, so clear it
        while(Serial1.read() != -1);
        // Write the buffer and wait for the CF to finish responding
        Serial1.write(req_buffer, 60);
        // Read into buffer and convert to a response
        if (Serial1.readBytes(res_buffer, 60) < 60) {
            last_connected ++;
            continue;
        }
        last_connected = 0;
        memcpy(&res, res_buffer, 60);
        // if (Serial) Serial.print(res.px);
        // if (Serial) Serial.print(res.py);
        // if (Serial) Serial.println(res.pz);
        // if (Serial) Serial.print(res.vx);
        // if (Serial) Serial.print(res.vy);
        // if (Serial) Serial.println(res.vz);
        // if (Serial) Serial.print(res.ax);
        // if (Serial) Serial.print(res.ay);
        // if (Serial) Serial.println(res.az);
        // if (Serial) Serial.println(res.r);
        // if (Serial) Serial.println(res.p);
        // if (Serial) Serial.println(res.y);
        // if (Serial) Serial.print(res.yr);
        // if (Serial) Serial.print(res.pr);
        // if (Serial) Serial.println(res.rr);
        // if (Serial) Serial.println("MPC");
        if (!response_valid(&res)) {
            if (Serial) Serial.println("Response invalid");
            continue;
        }
        tinyVector x0 = (tinyVector(10) << (tinytype)res.px,(tinytype)res.py,(tinytype)res.pz,(tinytype)res.vx,(tinytype)res.vy,(tinytype)res.vz, (tinytype)res.ax, (tinytype)res.ay, (tinytype)res.az, (tinytype)res.y).finished();
        solver->work->x.col(0) = x0;        
        
        for (int i = 0; i < solver->work->N; i ++) {
            if (counter < end) {

                double x_orig, y_orig;
                int sweep_phase = counter / sweepperiod; 

                // // Position (fig 8)
                // ref(0, i) = 0.5 * sin(2 * freq * ((counter + 10) * DT + DT * i));
                // ref(1, i) = 0.5 * sin(freq * ((counter + 10) * DT + DT * i));
                // ref(2, i) = 0.5;  // Constant Z

                // // Position (fig s)
                // ref(0, i) = 0.5 * sin(1 * freq * ((counter + 10) * DT + DT * i));  
                // ref(1, i) = 0.5 * sin(3 * freq * ((counter + 10) * DT + DT * i));  
                // ref(2, i) = 0.5;
                
                // const double S = 2.0;  // Speed-up factor (adjust as needed)

                // if ((S * ((counter + 10) * DT + i * DT)) < period / 2) {
                //     // Move along X-axis
                //     x_orig = 0.5 * (S * ((counter + 10) * DT + i * DT));
                //     y_orig = 0.0;
                // } 
                // else if((S * ((counter + 10) * DT + i * DT)) < period) {
                //     // Move along Y-axis
                //     x_orig = 0.5 * (period / 2);
                //     y_orig = 0.5 * (S * ((counter + 10) * DT + i * DT) - (period / 2));
                // }
                // else {
                //     // Move along Z-axis
                //     x_orig = 0.5 * (period / 2);
                //     y_orig = 0.5 * (period / 2);    
                // }

                // if (((counter + 10) * DT + i * DT) < period / 2) {
                //     // Move along X-axis
                //     x_orig = 0.5 * ((counter + 10) * DT + i * DT);
                //     y_orig = 0.0;
                // } else {
                //     // Move along Y-axis
                //     x_orig = 0.5 * (period / 2);
                //     y_orig = 0.5 * (((counter + 10) * DT + i * DT) - (period / 2));
                // }

                // // position (L shaped)
                // ref(0, i) = cos_theta * x_orig - sin_theta * y_orig; 
                // ref(1, i) = sin_theta * x_orig + cos_theta * y_orig; 
                // ref(2, i) = 0.5;

                // // position (circle)
                // ref(0, i) = 0.5 * sin(freq * ((counter + 10) * DT + DT * i));  
                // ref(1, i) = 0.5 * cos(freq * ((counter + 10) * DT + DT * i));  
                // ref(2, i) = 0.5;
                // ref(2, i) = 0.5 + 0.2 * sin(freq * ((counter + 10) * DT + DT * i));  // Slanted Z
                
                // position (sweep)
                ref(0, i) = x_increment * sweep_phase;
                ref(1, i) = -0.5 * cos(freq * ((counter + 10) * DT + DT * i)) + 0.5; // Sweeping in y
                ref(2, i) = 0.2;  // Constant Z

                // // position (helix)
                // ref(0, i) = 0.5 * sin(freq * ((counter + 10) * DT + DT * i));  
                // ref(1, i) = 0.5 * cos(freq * ((counter + 10) * DT + DT * i));  
                // ref(2, i) = 0.2 * ((counter + 10) * DT + DT * i); // Linear increase in height

                // if ((S * (counter * DT + i * DT)) < period / 2) {
                //     // Move along X-axis
                //     x_orig = 0.5 * (S * (counter * DT + i * DT));
                //     y_orig = 0.0;
                // } 
                // else if((S * (counter * DT + i * DT)) < period) {
                //     // Move along Y-axis
                //     x_orig = 0.5 * (period / 2);
                //     y_orig = 0.5 * (S * (counter * DT + i * DT) - (period / 2));
                // }
                // else {
                //     // Move along Z-axis
                //     x_orig = 0.5 * (period / 2);
                //     y_orig = 0.5 * (period / 2);
                // }

                // if ((counter * DT + i * DT) < period / 2) {
                //     // Move along X-axis
                //     x_orig = 0.5 * (counter * DT + i * DT);
                //     y_orig = 0.0;
                // } else {
                //     // Move along Y-axis
                //     x_orig = 0.5 * (period / 2);
                //     y_orig = 0.5 * ((counter * DT + i * DT) - (period / 2));
                // }

                // // Apply 2D Rotation (L Shape)
                // solver->work->Xref(0, i) = cos_theta * x_orig - sin_theta * y_orig;
                // solver->work->Xref(1, i) = sin_theta * x_orig + cos_theta * y_orig;
                // solver->work->Xref(2, i) = 0.5;  // Constant Z

                // // Position (circle)
                // solver->work->Xref(0, i) = 0.5 * sin(freq * ((counter) * DT + DT * i));
                // solver->work->Xref(1, i) = 0.5 * cos(freq * ((counter) * DT + DT * i)); 
                // solver->work->Xref(2, i) = 0.5;  // Constant Z
                // solver->work->Xref(2, i) = 0.5 + 0.2 * sin(freq * ((counter) * DT + DT * i));  // Slanted Z
                
                
                // position (sweep)
                solver->work->Xref(0, i) = x_increment * sweep_phase;
                solver->work->Xref(1, i) = -0.5 * cos(freq * ((counter) * DT + DT * i)) + 0.5; // Sweeping in y
                solver->work->Xref(2, i) = 0.2;  // Constant Z

                // // Position (fig 8)
                // solver->work->Xref(0, i) = 0.5 * sin(2 * freq * ((counter) * DT + DT * i));
                // solver->work->Xref(1, i) = 0.5 * sin(freq * ((counter) * DT + DT * i));
                // solver->work->Xref(2, i) = 0.5;  // Constant Z
                
                // // Position (fig s)
                // solver->work->Xref(0, i) = 0.5 * sin(1 * freq * ((counter) * DT + DT * i));  
                // solver->work->Xref(1, i) = 0.5 * sin(3 * freq * ((counter) * DT + DT * i));  
                // solver->work->Xref(2, i) = 0.5;

                // // position (helix)
                // solver->work->Xref(0, i) = 0.5 * sin(freq * ((counter) * DT + DT * i));  
                // solver->work->Xref(1, i) = 0.5 * cos(freq * ((counter) * DT + DT * i));  
                // solver->work->Xref(2, i) = 0.5 + 0.2 * ((counter) * DT + DT * i); // Linear increase in height
            }

            else {
                solver->work->Xref(3, i) = 0;
                solver->work->Xref(4, i) = 0;
                solver->work->Xref(5, i) = 0;

                solver->work->Xref(6, i) = 0;
                solver->work->Xref(7, i) = 0;
                solver->work->Xref(8, i) = 0;
            }
        }
        solve(solver);
         

        u_cmd(0) = solver->work->u.col(0)(0);
        u_cmd(1) = solver->work->u.col(0)(1);
        u_cmd(2) = solver->work->u.col(0)(2);
        u_cmd(0) = lpf2pApply(&cmdxfilt, u_cmd(0));
        u_cmd(1) = lpf2pApply(&cmdyfilt, u_cmd(1));
        u_cmd(2) = lpf2pApply(&cmdzfilt, u_cmd(2));
        a_cmd(0) = solver->work->x.col(1)(6);
        a_cmd(1) = solver->work->x.col(1)(7);
        a_cmd(2) = solver->work->x.col(1)(8);

        
        g << 0, 0, 9.81;
        ref_a << solver->work->x.col(1)(6), solver->work->x.col(1)(7), solver->work->x.col(1)(8);
        
        Eigen::VectorXd z_bar(solver->work->nx + 1);
        z_bar(0) = 1;
        z_bar.tail(solver->work->nx) = solver->work->x.col(1) - solver->work->z_star.col(1);
        // z_bar.tail(solver->work->nx) = x0 - solver->work->z_star.col(1);
        
        // drag_est << m_z_x[1].dot(z_bar), m_z_z[1].dot(z_bar);
        // drag_uncert << (L_hat_x[1] * z_bar).norm(), (L_hat_z[1] * z_bar).norm();
       

        // if (Serial) Serial.println("mx");
        // if (Serial) Serial.println(m_z_x[1].dot(z_bar));
        // if (Serial) Serial.println("my");
        // if (Serial) Serial.println(m_z_y[1].dot(z_bar));
        // if (Serial) Serial.println("mz");
        // if (Serial) Serial.println(m_z_z[1].dot(z_bar));

        // req.m9 = m_z_x[1].dot(z_bar);
        req.m9 = solver->work->Xref(2, 0);
        req.m10 = (L_hat_x[1] * z_bar).norm();
        req.m11 = m_z_z[1].dot(z_bar);
        req.m12 = (L_hat_z[1] * z_bar).norm();
        // req.m13 = m_z_z[1].dot(z_bar);
        // req.m14 = (L_hat_z[1] * z_bar).norm();
        req.m13 = x0(0);
        req.m14 = solver->work->Xref(0, 0);
        req.m7 = x0(1);
        req.m8 = solver->work->Xref(1, 0);
        req.m2 = x0(2);

        // req.m9 = x0(0);
        // req.m10 = solver->work->Xref(0, 0);
        // req.m11 = x0(1);
        // req.m12 = solver->work->Xref(1, 0);
        // req.m13 = x0(2);
        // req.m14 = solver->work->Xref(2, 0);
        
        if (solver->settings->en_state_bound == 1) {
            ref_a(0) = ref_a(0) + m_z_x[1].dot(z_bar);
            ref_a(1) = ref_a(1) + m_z_y[1].dot(z_bar);
            ref_a(2) = ref_a(2) + m_z_z[1].dot(z_bar);
        }
        target_accel = ref_a + g;
        // Eigen::Vector3d target_force = 0.036 * (ref_a + g);
        Eigen::Vector3d target_z_ax = target_accel.normalized();
        Eigen::Vector3d target_x_c(std::cos(0), std::sin(0), 0);
        Eigen::Vector3d target_y_ax = target_z_ax.cross(target_x_c).normalized();
        Eigen::Vector3d target_x_ax = target_y_ax.cross(target_z_ax);
        // Matrix3d target_rotation; 
        // target_rotation << target_x_ax, target_y_ax, target_z_ax;
        target_jerk << u_cmd(0), u_cmd(1), u_cmd(2);
        computeomega(&omega, target_x_ax, target_y_ax, target_accel, target_jerk);
        c << 0, 0, target_accel.norm();

        req.m0 = omega(0);
        req.m1 = -omega(1);
        // req.m2 = 0;
        req.m3 = ref_a(2);
        req.m4 = solver->work->x.col(1)(6);
        req.m5 = solver->work->x.col(1)(7);
        req.m6 = solver->work->x.col(1)(8);
        // req.m7 = solver->work->status;
        // req.m8 = solver->work->iter;
        
        counter++;
    }
    ::vTaskSuspend(nullptr);
}

void setup() {
    ::Serial.begin(115'200);
    solver = new TinySolver();
    solver->cache = new TinyCache();
    solver->settings = new TinySettings();
    solver->work = new TinyWorkspace();
    solver->settings->abs_pri_tol = (tinytype)0.0005; // Primal tolerance
    solver->settings->abs_dua_tol = (tinytype)0.0005; // Dual tolerance
    solver->settings->max_iter = iter_conf; // Max iterations
    solver->settings->check_termination = 1; // Check termination
    solver->settings->en_state_bound = 1; // Constrain state
    solver->settings->en_input_bound = 0; // Constrain input

    if (solver->settings->en_state_bound == 0){
        deleteFileIfExists("nstates_train_gp.csv");
        deleteFileIfExists("ndrags_train_gp.csv");
    }

    Q_glob << 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1000;


    solver->work = new TinyWorkspace();
    solver->work->nx = 10; // Number of states
    solver->work->nu = 4; // Number of inputs
    solver->work->N = 10; // Horizon length
    solver->work->x = tinyMatrix::Zero(solver->work->nx, solver->work->N);
    x_prev = tinyMatrix::Zero(solver->work->nx, solver->work->N);
    solver->work->u = tinyMatrix::Zero(solver->work->nu, solver->work->N-1);
    solver->work->d1 = tinyMatrix::Ones(1, solver->work->N);
    solver->work->d2 = tinyMatrix::Ones(1, solver->work->N);
    solver->work->d3 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->d4 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->d7 = tinyMatrix::Ones(1, solver->work->N);
    solver->work->d8 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->q = tinyMatrix::Zero(solver->work->nx, solver->work->N);
    solver->work->r = tinyMatrix::Zero(solver->work->nu, solver->work->N-1);
    solver->work->p = tinyMatrix::Zero(solver->work->nx, solver->work->N);
    solver->work->d = tinyMatrix::Zero(solver->work->nu, solver->work->N-1);
    solver->work->v = tinyMatrix::Zero(solver->work->nx, solver->work->N);
    solver->work->vnew = tinyMatrix::Zero(solver->work->nx, solver->work->N);
    solver->work->sd11 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd11_p = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd13 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd13_p = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd16 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd16_p = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd22 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd22_p = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd24 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd24_p = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd29 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd29_p = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd35 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd35_p = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd36 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd36_p = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd46 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd46_p = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd49 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd49_p = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd413 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd413_p = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd710 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd710_p = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd711 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd711_p = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd713 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd713_p = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd812 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd812_p = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd813 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sd813_p = tinyMatrix::Zero(1, solver->work->N);
    solver->work->sz1 = tinyMatrix::Zero(3, solver->work->N);
    solver->work->sz1_p = tinyMatrix::Zero(3, solver->work->N);
    solver->work->sz2 = tinyMatrix::Zero(3, solver->work->N);
    solver->work->sz2_p = tinyMatrix::Zero(3, solver->work->N);
    solver->work->sz3 = tinyMatrix::Zero(11, solver->work->N);
    solver->work->sz3_p = tinyMatrix::Zero(11, solver->work->N);
    solver->work->sz4 = tinyMatrix::Zero(11, solver->work->N);
    solver->work->sz4_p = tinyMatrix::Zero(11, solver->work->N);
    solver->work->sz5 = tinyMatrix::Zero(2, solver->work->N);
    solver->work->sz5_p = tinyMatrix::Zero(2, solver->work->N);
    solver->work->sz9 = tinyMatrix::Zero(solver->work->nx, solver->work->N);
    solver->work->sz9_p = tinyMatrix::Zero(solver->work->nx, solver->work->N);
    solver->work->sz10 = tinyMatrix::Zero(3, solver->work->N);
    solver->work->sz10_p = tinyMatrix::Zero(3, solver->work->N);
    solver->work->sz11 = tinyMatrix::Zero(11, solver->work->N);
    solver->work->sz11_p = tinyMatrix::Zero(11, solver->work->N);
    solver->work->sz12 = tinyMatrix::Zero(2, solver->work->N);
    solver->work->sz12_p = tinyMatrix::Zero(2, solver->work->N);
    solver->work->szt1 = tinyMatrix::Zero(2, solver->work->N);
    solver->work->szt1_p = tinyMatrix::Zero(2, solver->work->N);
    solver->work->szt2 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->szt2_p = tinyMatrix::Zero(1, solver->work->N);
    solver->work->z = tinyMatrix::Zero(solver->work->nu, solver->work->N-1);
    solver->work->znew = tinyMatrix::Zero(solver->work->nu, solver->work->N-1);
    solver->work->vcnew = tinyMatrix::Zero(solver->work->nx, solver->work->N);
    solver->work->zcnew = tinyMatrix::Zero(solver->work->nu, solver->work->N-1);
    solver->work->g = tinyMatrix::Zero(solver->work->nx, solver->work->N);
    solver->work->y = tinyMatrix::Zero(solver->work->nu, solver->work->N-1);
    solver->work->ld11 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->ld13 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->ld16 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->ld22 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->ld24 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->ld29 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->ld35 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->ld36 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->ld46 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->ld49 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->ld413 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->ld710 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->ld711 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->ld713 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->ld812 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->ld813 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->lz1 = tinyMatrix::Zero(3, solver->work->N);
    solver->work->lz2 = tinyMatrix::Zero(3, solver->work->N);
    solver->work->lz3 = tinyMatrix::Zero(11, solver->work->N);
    solver->work->lz4 = tinyMatrix::Zero(11, solver->work->N);
    solver->work->lz5 = tinyMatrix::Zero(2, solver->work->N);
    solver->work->lz9 = tinyMatrix::Zero(solver->work->nx, solver->work->N);
    solver->work->lz10 = tinyMatrix::Zero(3, solver->work->N);
    solver->work->lz11 = tinyMatrix::Zero(11, solver->work->N);
    solver->work->lz12 = tinyMatrix::Zero(2, solver->work->N);
    solver->work->lzt1 = tinyMatrix::Zero(2, solver->work->N);
    solver->work->lzt2 = tinyMatrix::Zero(1, solver->work->N);
    solver->work->z_star = tinyMatrix::Zero(solver->work->nx, solver->work->N);
    solver->work->m_hat_z = tinyMatrix::Zero(1, solver->work->N);
    solver->work->H_z = tinyMatrix::Zero(10, solver->work->N);
    solver->work->Adyn = tinyMatrix::Zero(10, 10);
    solver->work->Bdyn = tinyMatrix::Zero(10, 4);
    solver->work->R = tinyVector::Zero(4);
    for (int i = 0; i < 10; i++)
        solver->work->Q[i] = (tinyMatrix(10, 10)  << 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1000).finished();

    solver->work->Adyn  << 1, 0, 0, 0.01, 0, 0, 5e-05, 0, 0, 0, 0, 1, 0, 0, 0.01, 0, 0, 5e-05, 0, 0, 0, 0, 1, 0, 0, 0.01, 0, 0, 5e-05, 0, 0, 0, 0, 1, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
    solver->work->Bdyn  << 1.66667e-07, 0, 0, 0, 0, 1.66667e-07, 0, 0, 0, 0, 1.66667e-07, 0, 5e-05, 0, 0, 0, 0, 5e-05, 0, 0, 0, 0, 5e-05, 0, 0.01, 0, 0, 0, 0, 0.01, 0, 0, 0, 0, 0.01, 0, 0, 0, 0, 0.01;
    solver->cache = new TinyCache();
    // solver->cache->rho = (tinytype)4e-05;

    for(int i = 0; i < 10; i++){
        solver->cache->Kinf[i] = tinyMatrix::Zero(4, 10);
        solver->cache->Pinf[i] = tinyMatrix::Zero(10, 10);
        solver->cache->Quu_inv[i] = tinyMatrix::Zero(4, 4);
        solver->cache->AmBKt[i] = tinyMatrix::Zero(10, 10);  
    }
    solver->cache->APf = (tinyVector::Zero(10));
    solver->cache->BPf = (tinyVector::Zero(4));

    double penalty = 0.01;
    solver->cache->rhox = penalty;
    solver->cache->rhod1 = 0.1;
    solver->cache->rhod2 = 0.1;
    solver->cache->rhod3 = 1;
    solver->cache->rhod4 = 1;
    solver->cache->rhod7 = 0.1;
    solver->cache->rhod8 = 1;
    

    solver->work->R  << 0.06, 0.06, 1e-05, 0.1;


    solver->cache->rhou = 9e-4;
    solver->cache->rhox0 = 0.2;
    
    MatrixXd Ktp1 = MatrixXd::Zero(solver->work->nu, solver->work->nx);
    MatrixXd Ptp1 = solver->cache->rhox0 * MatrixXd::Ones(solver->work->nx, 1).array().matrix().asDiagonal();
    MatrixXd _Kinf = MatrixXd::Zero(solver->work->nu, solver->work->nx);
    MatrixXd _Pinf = MatrixXd::Zero(solver->work->nx, solver->work->nx);
    MatrixXd R1 = solver->work->R.array().matrix().asDiagonal();
    MatrixXd Q1 = MatrixXd::Zero(solver->work->nx, solver->work->nx);
    Q1 = Q_glob;

    for (int i = 0; i < 1000; i++)
    {
        _Kinf = (R1 + solver->work->Bdyn.transpose() * Ptp1 * solver->work->Bdyn).inverse() * solver->work->Bdyn.transpose() * Ptp1 * solver->work->Adyn;
        _Pinf = Q1 + solver->work->Adyn.transpose() * Ptp1 * (solver->work->Adyn - solver->work->Bdyn * _Kinf);
        // if Kinf converges, break
        if ((_Kinf - Ktp1).cwiseAbs().maxCoeff() < 1e-5)
        {
            // std::cout << "Kinf converged after " << i + 1 << " iterations" << std::endl;
            break;
        }
        Ktp1 = _Kinf;
        Ptp1 = _Pinf;
    }


    // Compute cached matrices
    MatrixXd _Quu_inv = (R1 + solver->work->Bdyn.transpose() * _Pinf * solver->work->Bdyn).inverse();
    MatrixXd _AmBKt = (solver->work->Adyn - solver->work->Bdyn * _Kinf).transpose();
    for (int i = 0; i < 10; i ++) {
        solver->work->Q[i] = Q1;
        solver->cache->Kinf[i] = _Kinf;
        solver->cache->Pinf[i] = _Pinf;
        solver->cache->Quu_inv[i] = _Quu_inv;
        solver->cache->AmBKt[i] = _AmBKt;
    }

    P_N = solver->cache->Pinf[0];
    solver->work->fdyn = (tinyVector::Zero(solver->work->nx));
    solver->work->x_min = -100000000000000000.0*tinyMatrix::Ones(solver->work->nx, solver->work->N);
    solver->work->x_max = 100000000000000000.0*tinyMatrix::Ones(solver->work->nx, solver->work->N);
    solver->work->u_min = -100000000000000000.0*tinyMatrix::Ones(solver->work->nu, solver->work->N-1);
    solver->work->u_max = 100000000000000000.0*tinyMatrix::Ones(solver->work->nu, solver->work->N-1);
    solver->work->numStateCones = 0;
    solver->work->numInputCones = 0;
    solver->work->Xref = tinyMatrix::Zero(solver->work->nx, solver->work->N);
    solver->work->Uref = tinyMatrix::Zero(solver->work->nu, solver->work->N-1);
    solver->work->Qu = tinyVector::Zero(solver->work->nu);
    solver->work->primal_residual_state = (tinytype)0.0;
    solver->work->primal_residual_input = (tinytype)0.0;
    solver->work->dual_residual_state = (tinytype)0.0;
    solver->work->dual_residual_input = (tinytype)0.0;
    solver->work->status = 0;
    solver->work->iter = 0;
    solver->solution = new TinySolution();
    solver->solution->iter = 0;
    solver->solution->solved = 0;
    solver->solution->x = tinyMatrix::Zero(solver->work->nx, solver->work->N);
    solver->solution->u = tinyMatrix::Zero(solver->work->nu, solver->work->N-1);
   

    H_xz(0, 6) = 1;
    H_xz(1, 8) = 1;
    m_xz(1, 0) = 9.81;
    H_x(0, 6) = 1;
    z_bar(0) = 1;
    z_bar.tail(solver->work->nx).setZero();
    
    // return 0;
    SOC c_xyz; // 0
    c_xyz.A = H_xyz;
    c_xyz.b = m_xyz * z_bar;
    solver->work->socx.push_back(c_xyz);

    SOC c_lx; // 1
    L_x_r = L_x.block(0, 1, 11, 10);
    c_lx.A = L_x_r;
    c_lx.b = L_x * z_bar;
    solver->work->socx.push_back(c_lx);

    SOC c_lz; // 2
    L_z_r = L_z.block(0, 1, 11, 10);
    c_lz.A = L_z_r;
    c_lz.b = L_z * z_bar;
    solver->work->socx.push_back(c_lz);

    H_xy.row(0) = H_xyz.row(0);
    H_xy.row(1) = H_xyz.row(1);
    m_xy.row(0) = m_xyz.row(0);
    m_xy.row(1) = m_xyz.row(1);

    SOC c_xy; // 3
    c_xy.A = H_xy;
    c_xy.b = m_xy * z_bar;
    solver->work->socx.push_back(c_xy);

    SOC c_ly; // 4
    L_y_r = L_y.block(0, 1, 11, 10);
    c_ly.A = L_y_r;
    c_ly.b = L_y * z_bar;
    solver->work->socx.push_back(c_ly);


    // newly added
    SOC soc1; // 5
    soc1.A = MatrixXd::Zero(2, 10);
    soc1.A(0, 6) = 1;
    soc1.A(1, 7) = 1;
    soc1.b = tinyVector::Zero(2);
    soc1.c = tinyVector::Zero(10);
    soc1.c(8) = tan(12 * M_PI / 180);
    soc1.d = tinyVector::Zero(1);
    soc1.d(0) = 9.81 * tan(12 * M_PI / 180);
    solver->work->socx.push_back(soc1);


    for (int j = 0; j < solver->work->N; j++){
         solver->work->SC.push_back(solver->work->socx);
    }    
    _p << 0, 0, 0;
    _v << 0, 0, 0;
    _a << 0, 0, 0;
    _j << 0, 0, 0;
    _w << 0, 0, 0; 
    c << 0, 0, 0; 
    _c << 0, 0, 0;
    tau << 0, 0, 0;
    ref_p << 0, 0, 0; 
    ref_v << 0, 0, 0; 
    Eigen::MatrixXd ref_R = Eigen::MatrixXd::Identity(3, 3);
    MatrixXd target_rotation = Eigen::MatrixXd::Identity(3, 3);

    P_COEFF_FOR << 10, 10, 10;
    I_COEFF_FOR << 0, 0, 0;
    D_COEFF_FOR << 10, 10, 15;

    P_COEFF_TOR << 0.1, 0.1, 0.1;
    I_COEFF_TOR << 0, 0, 0;
    D_COEFF_TOR << 0.01, 0.01, 0.01;


    P_COEFF_TOR << 1.0, 1.0, 1.0;
    I_COEFF_TOR << 0, 0, 0;
    D_COEFF_TOR << 0.1, 0.1, 0.1;


    for(int i = 0; i< 10; i++) {
        VectorXd a(11);
        MatrixXd b = MatrixXd::Zero(11, 11);
        m_z_x[i] = a;
        L_hat_x[i] = b;
        m_z_y[i] = a;
        L_hat_y[i] = b;
        m_z_z[i] = a;
        L_hat_z[i] = b;
    }

    req.m0 = 0.0f;
    req.m1 = 0.0f;
    req.m2 = 0.0f;
    req.m3 = 9.81f;
    req.m4 = 0.0f;
    req.m5 = 0.0f;
    req.m6 = 0.0f;
    req.m7 = 0.0f;
    req.m8 = 0.0f;
    req.m9 = 0.0f;
    req.m10 = 0.0f;
    req.m11 = 0.0f;
    req.m12 = 0.0f;
    req.m13 = 0.0f;
    req.m14 = 0.0f;
    lpf2pInit(&cmdxfilt, 100, 10);
    lpf2pInit(&cmdyfilt, 100, 10);
    lpf2pInit(&cmdzfilt, 100, 10);

    while (::millis() < 2'000) {
    }
    
    // ::xTaskCreate(DynamicsTask, "DynamicsTask", 1024, nullptr, 2, nullptr);
    ::xTaskCreate(GPTask, "GPTask", 2048, nullptr, 3, nullptr);
    ::xTaskCreate(ControllerTask, "ControllerTask", 4096, nullptr, 4, nullptr);

    ::Serial.println(PSTR("setup(): starting scheduler..."));
    ::Serial.flush();

    ::vTaskStartScheduler();
}

void loop() {}