#pragma once

// We need to define how the 2 talk over SPI.
// Teensy sends the motor commands and recieves the state

struct request { // 16
  float m0; // 4 bytes
  float m1; // 4 bytes
  float m2; // 4 bytes
  float m3; // 4 bytes
  float m4; // ax_cmd
  float m5; // ay_cmd
  float m6; // az_cmd
  float m7; // status
  float m8; // iter
  float m9; // mx
  float m10; // lx
  float m11; // my
  float m12; // ly
  float m13; // mz
  float m14; // lz
};

struct response { // 60 Bytes
  float px; // 4 bytes
  float py; // 4 bytes
  float pz; // 4 bytes
  float vx; // 4 bytes
  float vy; // 4 bytes
  float vz; // 4 bytes
  float ax; // 4 bytes
  float ay; // 4 bytes
  float az; // 4 bytes
  float y; // 4 bytes
  float p; // 4 bytes
  float r; // 4 bytes
  float yr; // 4 bytes
  float pr; // 4 bytes
  float rr; // 4 bytes
};
