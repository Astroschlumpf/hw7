#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

// ************* Funktionen fuer Ableitung F = {X', Y', Z'} ************* //
void vecDiv(double* F, double t, double* vec, double mu){
  // vec[] traegt {x, x', y, y', z, z'} in sich
  // F = {x', x'', y', y'', z', z''}
  // r = sqrt((vec[0] + mu)*(vec[0] + mu) + vec[2]*vec[2] + vec[4]*vec[4])
  // s = sqrt((vec[0] + mu - 1)*(vec[0] + mu - 1) + vec[2]*vec[2] + vec[4]*vec[4])
  
  F[0] = vec[1];
  F[1] = vec[0] + 2 * vec[3] - ((1 - mu) * (vec[0] + mu))/(((vec[0] + mu)*(vec[0] + mu) + vec[2]*vec[2] +
	 vec[4]*vec[4]) * sqrt((vec[0] + mu)*(vec[0] + mu) + vec[2]*vec[2] + vec[4]*vec[4])) -
	 (mu * (vec[0] - 1 + mu))/(((vec[0] + mu - 1)*(vec[0] + mu - 1) + vec[2]*vec[2] + vec[4]*vec[4]) * // worst ... equation ... ever ...
	 sqrt((vec[0] + mu - 1)*(vec[0] + mu - 1) + vec[2]*vec[2] + vec[4]*vec[4]));
  F[2] = vec[3];
  F[3] = vec[2] - 2 * vec[1] - ((1 - mu) * vec[2])/(((vec[0] + mu)*(vec[0] + mu) + vec[2]*vec[2] +
         vec[4]*vec[4]) * sqrt((vec[0] + mu)*(vec[0] + mu) + vec[2]*vec[2] + vec[4]*vec[4])) - (mu *
         vec[2])/(((vec[0] + mu - 1)*(vec[0] + mu - 1) + vec[2]*vec[2] + vec[4]*vec[4]) * sqrt((vec[0] + // oh wait... one more of these...
         mu - 1)*(vec[0] + mu - 1) + vec[2]*vec[2] + vec[4]*vec[4]));
  F[4] = vec[5];
  F[5] = -((1 - mu) * vec[4])/(((vec[0] + mu)*(vec[0] + mu) + vec[2]*vec[2] + vec[4]*vec[4]) *
         sqrt((vec[0] + mu)*(vec[0] + mu) + vec[2]*vec[2] + vec[4]*vec[4])) - (mu *
         vec[4])/(((vec[0] + mu - 1)*(vec[0] + mu - 1) + vec[2]*vec[2] + vec[4]*vec[4]) * // Immerhin keinen Speicherplatz fuer r und s blockiert!
         sqrt((vec[0] + mu - 1)*(vec[0] + mu - 1) + vec[2]*vec[2] + vec[4]*vec[4]));
}
// ********************************************************************* //

void mimaNorm(double* tvec, double* mvec){
  tvec[0] = max(abs(tvec[0]-mvec[0]),abs(tvec[1]-mvec[1]));
  tvec[0] = max(tvec[0],abs(tvec[2]-mvec[2]));
  tvec[0] = max(tvec[0],abs(tvec[3]-mvec[3]));
  tvec[0] = max(tvec[0],abs(tvec[4]-mvec[4]));
  // Max. war entweder in [0], [1], ..., oder wird nun durch [5] ueberschrieben
  tvec[0] = max(tvec[0],abs(tvec[5]-mvec[5]));
  // cout << tvec[0] << endl; // debug
}

// *************** Dormand-Prince, Runge-Kutta + step size ************* //
void vecDoPr(void(*vecDiv)(double*, double, double*, double),
	     double* dis, double t, double* vec, double mu, double epsi){
  //
  // btRK45[5][5] = {0,    |  0,            0,             0,            0,          0,              0,          0,
  //                 0.2,  |  0.2,          0,             0,            0,          0,              0,          0,
  //                 0.3,  |  0.075,        0.225,         0,            0,          0,              0,          0,
  //                 0.8,  |  44/45.,       -56/15.,       32/9.,        0,          0,              0,          0,
  //                 8/9., |  19372/6561.,  -25360/2187.,  64448/6561.,  -212/729.,  0,              0,          0,
  //                 1,    |  9017/3168.,   -355/33.,      46732/5247.,  49/176.,    -5103/18656.,   0,          0,
  //                 1,    |  35/384.,      0,             500/1113.,    125/192.,   -2187/6784.,    11/84.,     0,
  //                ------------------------------------------------------------------------------------------------
  //                 0,    |  35/384.,      0,            500/1113.,    125/192.,   -2187/6784.,     11/84.,     0,
  //                 0,    |  5179/57600.,  0,            7571/16695.,  393/640.,   -92097/339200.,  187/2100.,  0.025};
  //
  
  double kA[6] = {0.,0.,0.,0.,0.,0.}; // 3 Dimensionen trotz Nichtbeachten fuer z
  double kB[6] = {0.,0.,0.,0.,0.,0.}; // Kopie 1 der urspruenglichen Werte == NICHT ueberschreiben
  double kC[6] = {0.,0.,0.,0.,0.,0.}; // Kopie 2 der urspr...
  double kD[6] = {0.,0.,0.,0.,0.,0.}; // Kopie 3 der urspr...
  double kE[6] = {0.,0.,0.,0.,0.,0.}; // Kopie 4 der urspr...
  double kF[6] = {0.,0.,0.,0.,0.,0.}; // Kopie 5 der urspr...
  double kG[6] = {0.,0.,0.,0.,0.,0.}; // Kopie 6 der urspr... -> "RK7"
  double tmpV[6]; // Kopie zu "RK6" fuer anschliessende dt-Berechnung
  // Die Konstruktion ueber sechs verschiedene k[] kann sicher optimiert werden,
  // mir fehlte bei jedem meiner Versuche dabei die Uebersicht

  for(int i = 0; i < 6; i++){
    kA[i] = vec[i]; // initiale Werte in kA
    tmpV[i] = vec[i]; // Kopie der Yps(n) fuer Schritt Yps(n+1) als "RK6"-Berechnung
  }
  vecDiv(kA, t, vec, mu); // k1
  
  for(int i = 0; i < 6; i++){
    kB[i] = vec[i] + 0.2 * dis[0] * kA[i];
    tmpV[i] += dis[0] * (35/384.) * kA[i]; // k1 hinzufuegen
    vec[i] += dis[0] * (5179/57600.) * kA[i];
  }
  vecDiv(kB, t+0.2*dis[0], kB, mu); // k2: t+0.2*dt trotz Zeitunabhaengigkeit

  for(int i = 0; i < 6; i++){
    kC[i] = vec[i] + dis[0] * (0.075 * kA[i] + 0.225  * kB[i]);
    // vec[i]
    // tmpV[i] kB[] == k2 fliesst weder in Yps6 noch Yps7 ein
  }
  vecDiv(kC, t+0.3*dis[0], kC, mu); // k3
  
  for(int i = 0; i < 6; i++){
    kD[i] = vec[i] + dis[0] * ((44/45.) * kA[i] - (56/15.) * kB[i] + (32/9.) * kC[i]);
    tmpV[i] += dis[0] * (500/1113.) * kC[i]; // k3 hinzufuegen
    vec[i] += dis[0] * (7571/16695.) * kC[i];
  }
  vecDiv(kD, t+0.8*dis[0], kD, mu); // k4
  
  for(int i = 0; i < 6; i++){
    kE[i] = vec[i] + dis[0] * ((19372/6561.) * kA[i] - (25360/2187.) * kB[i] + (64448/6561.) * kC[i] - (212/729.) * kD[i]);
    tmpV[i] += dis[0] * (125/192.) * kD[i]; // k4 hinzufuegen
    vec[i] += dis[0] * (393/640.) * kD[i];
  }
  vecDiv(kE, t+(8/9.)*dis[0], kE, mu); // k5
  
  for(int i = 0; i < 6; i++){
    kF[i] = vec[i] + dis[0] * ((9017/3168.) * kA[i] - (355/33.) * kB[i] + (46732/5247.) * kC[i] + (49/176.) * kD[i] - (5103/18656.) * kE[i]);
    tmpV[i] -= dis[0] * (2187/6784.) * kE[i]; // k5 hinzufuegen
    vec[i] -= dis[0] * (92097/339200.) * kE[i]; // Vorzeichen beachten!
  }
  vecDiv(kF, t+dis[0], kF, mu); // k6
  
  for(int i = 0; i < 6; i++){
    kG[i] = vec[i] + dis[0] * ((35/384.) * kA[i] + (500/1113.) * kC[i] + (125/192.) * kD[i] - (2187/6784.) * kE[i] + (11/84.) * kF[i]); // + 0*kB[] !!
    tmpV[i] += dis[0] * (11/84.) * kF[i]; // k6 hinzufuegen
    vec[i] += dis[0] * (187/2100.) * kF[i];
  }
  vecDiv(kG, t+dis[0], kG, mu); // k7

  for(int i = 0; i < 6; i++){
    vec[i] += dis[0] * 0.025 * kG[i]; // kG[] == k7 hinzufuegen
    // kein tmpV[] mehr, enthaelt Yps6
  }

  // cout << tmpV[0]-vec[0] << endl; // debug

  mimaNorm(tmpV,vec); // Maximum der Differenzen (tmpV-vec) in tmpV[0] schreiben
  dis[0] = dis[0]  * pow((epsi/tmpV[0]), 0.2); // neues dt, in tmpV[0] wurde zuvor in mimaNorm     * 0.8
                                               // das Maximum aller Differenzen gespeichert
                                               // und in die referenzierte Variable dis == dt geschrieben

  // ANMERKUNG: an dieser Stelle wird dt EXTREM klein, aehnlich dt -> 0 fuer n -> oo
  // Des Weiteren ergibt tmpV[0] - vec[0] einen Wert von etwa 1.2 und wird in keinem
  // Schritt verkleinert. Die Funktion mimaNorm(d*,d*) scheint aber korrekt zu
  // arbeiten, was ich in einem externen Programm getestet habe. In Kombination
  // fuehren die (fehlerhaften) Berechnungen dazu, dass meine Ausgabe nur asymptotische
  // Werte in x und y ausgibt, was nicht die gesuchte Bahn ist. Ich konnte aber keinen
  // Fehler im Code ausmachen (ich uebersehe ihn also immer).
}
// ********************************************************************* //

int main(void){
  const double MU = 7.349/597.4;
  double t = 0.0; // Start: 0.0, Ende: 35 (knapp ueber 2 Perioden)
  int ind = 0; // int N = 1;

  double ord;
  cout << "Fehlertoleranz-Groessenordnung (3 bis 10):" << endl;
  cin >> ord;
  const double EPS = pow(10, -ord);

  double dt = 0.001; // initiale Schrittlaenge -> spaeter dt = dt * pow((EPS/|(y_6 - y_7)|), 0.2) (nach Vorlesung)

  // cout << EPS << endl; // debug
  
  ofstream Ausg("arenstorf.csv");

  double Yps[6]; // Yps = {[x,x'][y,y'][z,z']}
  // initiale Werte x(0), y(0) und z(0)
  Yps[0] = 0.994; Yps[1] = 0.; Yps[2] = 0.; Yps[3] = -2.00158510637908; Yps[4] = 0.; Yps[5] = 0.;
  
  Ausg << t << "\t" << Yps[0] << "\t" << Yps[2] << endl; // t, x, y; ohne z

  while((t < 17.1) && (ind < 3000000)){ // 35 als Endpunkt, ind als Fehlerkorrektur
    vecDoPr(vecDiv, &dt, t, Yps, MU, EPS); // Runge-Kutta aufrufen inklusive Schrittlaengenberechnung   
    // cout << ind << endl; // debug
    // cout << t << "\t" << Yps[0] << "\t" << Yps[2] << endl; // t | x | y   ausgeben
    Ausg << t << "\t" << Yps[0] << "\t" << Yps[2] << endl; // t | x | y   ausgeben
    t += dt;
    ind++;
  }
  
  cout << ind << " Schritte geschrieben." << endl; // debug
         
  Ausg.close();
  return 0;
}