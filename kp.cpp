#include <iostream>
#include <cmath>
#include "little.h"
#define M_PI 3.14159265358979323846
int main()
{
    int h = 154 * 1e3;
    double hag = 0.1;
    double P = 10.25 * 1e3;
    double Wisl = 3355;
    double betta = P / Wisl;
    double dMtmax = 1500;
    double m0 = 3015;
    double tv = 14;
    double tmax = dMtmax / betta;
    double t1 = 0.8 * tmax;
    double dtetta1 = -M_PI / 2 / (tmax - tv);
    double tetta2 = -M_PI / 6;
    double t2 = t1 + 50;
    rungekutta(t1, t2, dtetta1, tetta2, m0, h, P, Wisl, hag);
    control pop = kosnet(t1, t2, dtetta1, tetta2, m0, h, P, Wisl, hag);
    rungekutta(t1, t2, pop.dtetta1, pop.tetta2, m0, h, P, Wisl, hag, 1);
    std::cout << pop.dtetta1 << '\t' << pop.tetta2 << '\n';
    control flip = vremya(pop, pop.dtetta1, pop.tetta2, m0, h, P, Wisl, hag, 1);
    std::cout << flip.t1 << '\t' << flip.t2 << '\t' << flip.dtetta1 << '\t' << flip.tetta2 << '\n';
}
#include <iostream>
#include <cmath>
#include "little.h"
#include "C:\Users\trot5\source\repos\Интегрирование орбиты\Интегрирование орбиты\kepler.h"
#include <fstream>
#include <vector>
#include <iomanip>
#include <string>
#define M_PI 3.14159265358979323846
// P = 10.25 для 1 траектории h=154 lymbda 2.5 для 2 траектории h=307 lymbda 2.5
// P = 8.55 для 1 траектории h=154 lymbda 2.75 для 2 траектории h=198 lymbda 3.75
using namespace std;
int R = 1738 * 1e3;       // Радиус Луны (м)
double Mu = 4.903 * 1e12; // Гравитационная постоянная (км^3/с^2)
double dMtmax = 1500;     // (кг) Предельная масса рабочего запаса топлива ДУ СВ
double Mcons = 615;       // (кг) Масса конструкции СВ
double tv = 14;           // (с) Продолжительность вертикального участка выведения
double nachhag = 0.1;

double dVxdt(double t, double x, double y, double t1, double t2, double dtetta1, double tetta2, double m0, double P, double betta, bool j) // Производная Скорости по координате х
{
    return fP(t, t1, t2, P, j) * cos(tangag(t, t1, t2, dtetta1, tetta2, j)) / dmdt(m0, betta, t, t1, t2, j) + gx(x, y);
}
double dVydt(double t, double x, double y, double t1, double t2, double dtetta1, double tetta2, double m0, double P, double betta, bool j) // Производная скорости по координате y
{
    return fP(t, t1, t2, P, j) * sin(tangag(t, t1, t2, dtetta1, tetta2, j)) / dmdt(m0, betta, t, t1, t2, j) + gy(x, y);
}
double fP(double t, double t1, double t2, double P, bool j)
{
    if (j == true)
    {
        return P;
    }
    else
    {
        return 0;
    }
}
double dxdt(double Vx) // Произовдная по координате х
{
    return Vx;
}
double dydt(double Vy) // Производная по координате у
{
    return Vy;
}
double dmdt(double m0, double betta, double t, double t1, double t2, bool j) //  топливo
{
    if ((j == true) && (t <= t1))
    {
        return m0 - t * betta;
    }
    if (j == false)
    {
        return m0 - t1 * betta;
    }
    if ((t >= t2) && (j == true))
    {
        return m0 - t1 * betta - betta * (t - t2);
    }
}
double gx(double x, double y)
{
    return -(Mu * x) / pow((x * x + pow((R + y), 2)), 1.5);
}
double gy(double x, double y)
{
    return -(Mu * (R + y)) / pow((x * x + pow((R + y), 2)), 1.5);
}

double tangag(double t, double t1, double t2, double dtetta1, double tetta2, bool j) // Тангаж
{
    if ((0 <= t) && (t <= tv))
    {
        return M_PI / 2;
    }
    if ((j == true) && (t < t2 - 0.1))
    {
        return M_PI / 2 + dtetta1 * (t - tv);
    }
    if (j == false)
    {
        return M_PI / 2 + dtetta1 * (t1 - tv);
    }
    if ((t >= t2) && (j == true))
    {
        return tetta2;
    }
}
double uglnaklonakstarkgor(double Vy, double Vx) // угол наклона к стартовому горизонту
{
    if (Vx == 0)
        return M_PI / 2;
    else
        return atan2(Vy, Vx);
}
double fi(double Tetta, double uglnaklona)
{
    return Tetta - uglnaklona;
}
double Risl(double h)
{
    return R + h;
}
double Visl(double h)
{
    double f = sqrt(Mu / Risl(h));
    return sqrt(Mu / Risl(h));
}

double scalyartimes(double x, double y, double Vx, double Vy)
{
    return (x * Vx + (y + R) * Vy);
}
double sTetta(double x, double y, double Vx, double Vy, double r, double V) // Угол наклона к местному горизонту
{
    double k = scalyartimes(x, y, Vx, Vy);
    if (r == R)
        return M_PI / 2;
    else
        return asin(k / (r * V));
}
double rras(double x, double y)
{
    return sqrt(pow(x, 2) + pow(y + R, 2));
}
bool flagept(double t, double t1, double t2)
{
    if (t < t1)
    {
        return 1;
    }
    if ((t >= t1) && (t < t2))
    {
        return 0;
    }
    if (t >= t2)
    {
        return 1;
    }
}
double dVypravdt(double t, double t1, double t2, double m, double P, double alf)
{
    return P / m * (1 - cos(alf));
}
double time(double t, double hag, double t1, double t2, double i, double V, double V_orb)
{
    if (((t + hag) > t1) && (abs(t - t1) < hag))
    {
        double dt1 = t1 - t;
        if (t - t1 < 0)
        {
            return dt1;
        }
        if (abs(t - t1) < 0.001)
        {
            dt1 = t1 - hag + nachhag - t1;
            double dt2 = t - dt1 + hag - t1;
            return dt1;
        }
        if (abs(i + hag - t) < 0.001)
        {
            return hag;
        }
    }
    if (((t + hag) > t2) && (abs(t - t2) < hag))
    {
        double dt1 = t2 - t;
        if (t - t2 < 0)
        {
            return dt1;
        }
        if (abs(t - t2) < 0.001)
        {
            dt1 = t2 - hag + nachhag - t2;
            double dt2 = t - dt1 + hag - t2;
            return dt1;
        }
        if (abs(i + hag - t) < 0.001)
        {
            return hag;
        }
    }
    if (abs(V - V_orb) < 2)
        return hag;
    if (V < V_orb)
    {
        return nachhag;
    }
}
double dVgdt(double x, double y, double Vx, double Vy)
{
    double g = sqrt(pow(gx(x, y), 2) + pow(gy(x, y), 2));
    double r = rras(x, y);
    double V = sqrt(Vx * Vx + Vy * Vy);
    double Tetta = M_PI / 2 - sTetta(x, y, Vx, Vy, r, V);
    return g * cos(Tetta);
}
lol rungekutta(double t1, double t2, double dtetta1, double tetta2, double m0, double h, double P, double Wisl, double hag, bool isNeedtoPrint)
{
    kepler okm;
    fstream pleg, plug;
    if (isNeedtoPrint)
    {
        string asd1, asd2;
        cout << "Vvedite imya faila" << '\n';
        cin >> asd1;
        cout << "Vvedite imya vtorogo faila" << '\n';
        cin >> asd2;
        asd1 += ".csv";
        asd2 += ".csv";
        pleg.open(asd1, ios_base::out);
        plug.open(asd2, ios_base::out);
    }
    if (pleg.fail())
    {
        cout << "idiot" << '\n';
        exit(1);
    }
    double t = 0;
    double tk = 0;
    double V_orb = Visl(h);
    double betta = P / Wisl;
    if (isNeedtoPrint)
    {
        pleg << "Производная по массе" << '\n';
        pleg << proizvm(t1, t2, dtetta1, tetta2, m0, h, P, Wisl, hag) << '\n';
        pleg << "Производная по скорости истечения" << '\n';
        pleg << proizvWisl(t1, t2, dtetta1, tetta2, m0, h, P, Wisl, hag) << '\n';
        pleg << "Управляющие параметры" << '\n';
        pleg << "t1" << ';' << "t2" << ';' << "dtetta1" << ';' << "tetta2" << '\n';
        pleg << t1 << ';' << t2 << ';' << dtetta1 << ';' << tetta2 << '\n';
        pleg << "t" << ';' << "x" << ';' << "y" << ';' << "Vx" << ';' << "Vy" << ';' << "m" << ';' << "r" << ';' << "V" << ';' << "h" << ';' << "Tetta" << ';' << "Тангаж" << ';' << "Угол наклона к стартовому горизонту" << ';' << "alfa" << ';' << "fi" << '\n';
        plug << "t" << ';' << "x" << ';' << "y" << ';' << "Vx" << ';' << "Vy" << ';' << "m" << ';' << "r" << ';' << "V" << ';' << "h" << ';' << "Tetta" << ';' << "Тангаж" << ';' << "Угол наклона к стартовому горизонту" << ';' << "alfa" << ';' << "fi" << '\n';
    }
    double x_old = 0;
    double y_old = 0;
    double Vx_old = 0;
    double Vy_old = 0;
    double m_old = m0;
    double V = 0;
    double x_new = 0;
    double y_new = 0;
    double Vx_new = 0;
    double Vy_new = 0;
    double Vyprav_old = 0;
    double Vyprav1 = 0;
    double Vyprav2 = 0;
    double Vyprav3 = 0;
    double Vg_old = 0;
    double Vg1 = 0;
    double Vyprav_new = 0;
    double Vts1 = 0;
    double Vts2 = 0;
    double x;
    double y;
    double Vx;
    double Vy;
    double r = rras(x_old, y_old);
    double g = gx(x_old, y_old);
    double Vg2 = 0;
    double Vg3 = 0;
    double Vg_new = 0;
    double m1aut = 0;
    double m2put = 0;
    double Vts3 = 0;
    int a = 0;
    lol poh;
    poh.r = 0;
    poh.stetta = 0;
    poh.V = 0;
    poh.m = 0;
    bool j = 1;
    double i = 0;
    double Tetta = sTetta(x_old, y_old, Vx_old, Vy_old, r, V);
    double tang = tangag(t, t1, t2, dtetta1, tetta2, j);
    double uglnaklona = uglnaklonakstarkgor(Vy_old, Vx_old);
    vector<double> ti;
    if (isNeedtoPrint)
        pleg << t << ';' << x_old << ';' << y_old << ';' << Vx_old << ';' << Vy_old << ';' << m0 << ';' << r << ';' << V << ';' << r - R << ';' << 180 / M_PI * Tetta << ';' << 180 / M_PI * tang << ';' << 180 / M_PI * uglnaklona << ';' << 180 / M_PI * (tang - uglnaklona) << ';' << fi(Tetta, uglnaklona) * 180 / M_PI << '\n';
    plug << t << ';' << x_old << ';' << y_old << ';' << Vx_old << ';' << Vy_old << ';' << m0 << ';' << r << ';' << V << ';' << r - R << ';' << 180 / M_PI * Tetta << ';' << 180 / M_PI * tang << ';' << 180 / M_PI * uglnaklona << ';' << 180 / M_PI * (tang - uglnaklona) << ';' << fi(Tetta, uglnaklona) * 180 / M_PI << '\n';
    while (true)
    {

        a += 1;
        if ((abs(t - t1) < hag) || (abs(t - t2) < hag))
        {
            ti.push_back(t);
            if (((t + hag) > t1) && (abs(t - t1) < hag))
            {
                i = ti.at(0);
                a -= 1;
            }
            else
            {
                a -= 1;
                i = ti.at(2);
            }
        }
        if (V < V_orb)
        {
            hag = time(t, hag, t1, t2, i, V, V_orb);
        }
        j = flagept(t, t1, t2);
        double k1x = hag * dxdt(Vx_old);
        double k1y = hag * dydt(Vy_old);
        double k1Vx = hag * dVxdt(t, x_old, y_old, t1, t2, dtetta1, tetta2, m0, P, betta, j);
        double k1Vy = hag * dVydt(t, x_old, y_old, t1, t2, dtetta1, tetta2, m0, P, betta, j);

        double k2x = hag * dxdt(Vx_old + k1Vx / 2);
        double k2y = hag * dydt(Vy_old + k1Vy / 2);
        double k2Vx = hag * dVxdt(t + hag / 2, x_old + k1x / 2, y_old + k1y / 2, t1, t2, dtetta1, tetta2, m0, P, betta, j);
        double k2Vy = hag * dVydt(t + hag / 2, x_old + k1x / 2, y_old + k1y / 2, t1, t2, dtetta1, tetta2, m0, P, betta, j);

        double k3x = hag * dxdt(Vx_old + k2Vx / 2);
        double k3y = hag * dydt(Vy_old + k2Vy / 2);
        double k3Vx = hag * dVxdt(t + hag / 2, x_old + k2x / 2, y_old + k2y / 2, t1, t2, dtetta1, tetta2, m0, P, betta, j);
        double k3Vy = hag * dVydt(t + hag / 2, x_old + k2x / 2, y_old + k2y / 2, t1, t2, dtetta1, tetta2, m0, P, betta, j);

        double k4x = hag * dxdt(Vx_old + k3Vx);
        double k4y = hag * dydt(Vy_old + k3Vy);
        double k4Vx = hag * dVxdt(t + hag, x_old + k3x, y_old + k3y, t1, t2, dtetta1, tetta2, m0, P, betta, j);
        double k4Vy = hag * dVydt(t + hag, x_old + k3x, y_old + k3y, t1, t2, dtetta1, tetta2, m0, P, betta, j);

        double x_new = x_old + (k1x + 2 * k2x + 2 * k3x + k4x) / 6;
        double y_new = y_old + (k1y + 2 * k2y + 2 * k3y + k4y) / 6;
        double Vx_new = Vx_old + (k1Vx + 2 * k2Vx + 2 * k3Vx + k4Vx) / 6;
        double Vy_new = Vy_old + (k1Vy + 2 * k2Vy + 2 * k3Vy + k4Vy) / 6;

        if (isNeedtoPrint)
        {
            double k1yprav = hag * dVypravdt(t, t1, t2, dmdt(m0, betta, t, t1, t2, j), fP(t, t1, t2, P, j), tangag(t, t1, t2, dtetta1, tetta2, j) - uglnaklonakstarkgor(Vy_old, Vx_old));
            double k1g = hag * dVgdt(x_old, y_old, Vx_old, Vy_old);
            double k2yprav = hag * dVypravdt(t + hag / 2, t1, t2, dmdt(m0, betta, t + hag / 2, t1, t2, j), fP(t + hag / 2, t1, t2, P, j), tangag(t + hag / 2, t1, t2, dtetta1, tetta2, j) - uglnaklonakstarkgor(Vy_old + k1Vy / 2, Vx_old + k1Vx / 2));
            double k2g = hag * dVgdt(x_old + k1x / 2, y_old + k1y / 2, Vx_old + k1Vx / 2, Vy_old + k1Vy / 2);
            double k3yprav = hag * dVypravdt(t + hag / 2, t1, t2, dmdt(m0, betta, t + hag / 2, t1, t2, j), fP(t + hag / 2, t1, t2, P, j), tangag(t + hag / 2, t1, t2, dtetta1, tetta2, j) - uglnaklonakstarkgor(Vy_old + k2Vy / 2, Vx_old + k2Vx / 2));
            double k3g = hag * dVgdt(x_old + k2x / 2, y_old + k2y / 2, Vx_old + k2Vx / 2, Vy_old + k2Vy / 2);
            double k4yprav = hag * dVypravdt(t + hag, t1, t2, dmdt(m0, betta, t + hag, t1, t2, j), fP(t + hag, t1, t2, P, j), tangag(t + hag, t1, t2, dtetta1, tetta2, j) - uglnaklonakstarkgor(Vy_old + k3Vy, Vx_old + k3Vx));
            double k4g = hag * dVgdt(x_old + k3x, y_old + k3y, Vx_old + k3Vx, Vy_old + k3Vy);
            Vyprav_new = Vyprav_old + (k1yprav + 2 * k2yprav + 2 * k3yprav + k4yprav) / 6;
            Vg_new = Vg_old + (k1g + 2 * k2g + 2 * k3g + k4g) / 6;
        }
        V = sqrt(pow(Vx_new, 2) + pow(Vy_new, 2));
        if (V > V_orb)
        {
            hag /= 2;
            a = 0;
            continue;
        }
        t += hag;
        if ((abs(t - t1 - hag) < 0.0001) || (abs(t - t2 - hag) < 0.0001))
        {
            a += 1;
        }
        g = sqrt(pow(gx(x_new, y_new), 2) + pow(gy(x_new, y_new), 2));
        double m = dmdt(m0, betta, t, t1, t2, j);
        x_old = x_new;
        y_old = y_new;
        Vx_old = Vx_new;
        Vy_old = Vy_new;
        Vyprav_old = Vyprav_new;
        Vg_old = Vg_new;
        if (isNeedtoPrint)
        {
            if (abs(t - t1) < 0.001)
            {
                Vyprav1 = Vyprav_new;
                Vg1 = Vg_new;
                m1aut = m;
                Vts1 = Wisl * log(m0 / m) - Vg1 - Vyprav1;
            }
            if (abs(t - t2) < 0.001)
            {
                Vyprav2 = Vyprav_new;
                Vg2 = Vg_new;
                Vts2 = Vts1 - (Vg2 - Vg1) - (Vyprav2 - Vyprav1);
            }
            if (abs(V_orb - V) < 1e-6)
            {
                Vyprav3 = Vyprav_new;
                Vg3 = Vg_new;
                Vts3 = Vts2 + Wisl * log(m1aut / m) - (Vyprav3 - Vyprav2) - (Vg3 - Vg2);
            }
        }
        x = x_new;
        y = y_new;
        Vx = Vx_new;
        Vy = Vy_new;
        r = rras(x, y);
        double stetta = sTetta(x, y, Vx, Vy, r, V);
        poh.r = r;
        poh.V = V;
        poh.stetta = stetta;
        poh.m = m;
        Tetta = sTetta(x_new, y_new, Vx_new, Vy_new, r, V);
        tang = tangag(t, t1, t2, dtetta1, tetta2, j);
        uglnaklona = uglnaklonakstarkgor(Vy_new, Vx_new);
        if ((isNeedtoPrint))
        {
            plug << setprecision(15) << t << ';' << x_new << ';' << y_new << ';' << Vx_new << ';' << Vy_new << ';' << poh.m << ';' << r << ';' << V << ';' << r - R << ';' << 180 / M_PI * Tetta << ';' << 180 / M_PI * tang << ';' << 180 / M_PI * uglnaklona << ';' << 180 / M_PI * (tang - uglnaklona) << ';' << fi(Tetta, uglnaklona) * 180 / M_PI << ';' << g << '\n';
            if ((abs(V_orb - V) < 1e-6) || (a % 500 == 0) || (abs(t - tv) < 1e-3) || (abs(t - t1) < 1e-3) || (abs(t - t2) < 1e-3))
                pleg << setprecision(15) << t << ';' << x_new << ';' << y_new << ';' << Vx_new << ';' << Vy_new << ';' << poh.m << ';' << r << ';' << V << ';' << r - R << ';' << 180 / M_PI * Tetta << ';' << 180 / M_PI * tang << ';' << 180 / M_PI * uglnaklona << ';' << 180 / M_PI * (tang - uglnaklona) << ';' << fi(Tetta, uglnaklona) * 180 / M_PI << ';' << g << '\n';
        }
        if ((abs(t - t1) < 1e-3))
        {
            okm = dectokepler(x_new / 1000, (y_new + R) / 1000, 0, Vx_new / 1000, Vy_new / 1000, 0);
        }
        if (abs(V_orb - V) < 1e-6)
        {
            if (isNeedtoPrint)
            {
                pleg << "Потери скорости на управление" << '\n';
                pleg << "1-АУТ" << ';' << "ПУТ" << ';' << "2-АУТ" << '\n';
                pleg << Vyprav1 << ';' << Vyprav2 - Vyprav1 << ';' << Vyprav3 - Vyprav2 << '\n';
                pleg << "Гравитационные потери скорости" << '\n';
                pleg << "1-АУТ" << ';' << "ПУТ" << ';' << "2-АУТ" << '\n';
                pleg << Vg1 << ';' << Vg2 - Vg1 << ';' << Vg3 - Vg2 << '\n';
                pleg << "Характерестическая скорость" << '\n';
                pleg << Vts1 << ';' << Vts2 << ';' << Vts3 << '\n';
                pop znach = Hpa(okm.a, okm.e);
                pleg << "Параметры обриты" << '\n';
                pleg << "a" << ';' << "e" << ';' << "Ra" << ';' << "Rp" << '\n';
                pleg << okm.a << ';' << okm.e << ';' << znach.Ha << ';' << znach.Hp << '\n';
                pleg.close();
            }
            return poh;
        }
    }
}
control kosnet(double t1, double t2, double dtetta1, double tetta2, double m0, double h, double P, double Wisl, double hag)
{
    lol TT, yatoro;
    control Mira;
    Mira.dtetta1 = dtetta1;
    Mira.t1 = t1;
    Mira.t2 = t2;
    Mira.tetta2 = tetta2;
    TT = rungekutta(t1, t2, dtetta1, tetta2, m0, h, P, Wisl, hag);
    double R = Risl(h);
    double dr = TT.r - R;
    double dT = TT.stetta;
    double ddtetta1 = 1e-9;
    double dtetta2 = -1e-9;
    while (abs(dr) > 1e-2 || abs(dT) > (1e-4 * M_PI / 180))
    {
        lol TT1, TT2;
        TT1 = rungekutta(t1, t2, dtetta1 + ddtetta1, tetta2, m0, h, P, Wisl, hag);
        TT2 = rungekutta(t1, t2, dtetta1, tetta2 + dtetta2, m0, h, P, Wisl, hag);
        double dr1 = TT1.r - R;
        double dT1 = TT1.stetta;
        double ddrddt1 = (dr1 - dr) / ddtetta1;
        double ddTddt1 = (dT1 - dT) / ddtetta1;
        double dr2 = TT2.r - R;
        double dT2 = TT2.stetta;
        double ddrdt2 = (dr2 - dr) / dtetta2;
        double ddTdt2 = (dT2 - dT) / dtetta2;
        double det = ddrddt1 * ddTdt2 - ddrdt2 * ddTddt1;
        dtetta1 -= (ddTdt2 * dr - ddrdt2 * dT) / det;
        tetta2 -= (-ddTddt1 * dr + ddrddt1 * dT) / det;
        TT = rungekutta(t1, t2, dtetta1, tetta2, m0, h, P, Wisl, hag);
        dr = TT.r - R;
        dT = TT.stetta;
        // cout << dr << '\t' << dT << '\n';
        Mira.t1 = t1;
        Mira.t2 = t2;
        Mira.dtetta1 = dtetta1;
        Mira.tetta2 = tetta2;
    }
    return Mira;
}
double dmdt1(double t1, double t2, double dtetta1, double tetta2, double m0, double m, double h, double P, double Wisl, double hag, double eps)
{
    control Mira = kosnet(t1 + eps, t2, dtetta1, tetta2, m0, h, P, Wisl, hag);
    lol TT = rungekutta(t1 + eps, t2, Mira.dtetta1, Mira.tetta2, m0, h, P, Wisl, hag);
    return ((m0 - TT.m) - (m0 - m)) / (eps);
}
double dmdt2(double t1, double t2, double dtetta1, double tetta2, double m0, double m, double h, double P, double Wisl, double hag, double eps)
{
    control Mira = kosnet(t1, t2 + eps, dtetta1, tetta2, m0, h, P, Wisl, hag);
    lol TT = rungekutta(t1, t2 + eps, Mira.dtetta1, Mira.tetta2, m0, h, P, Wisl, hag);
    return ((m0 - TT.m) - (m0 - m)) / (eps);
}
control vremya(control Mira, double dtetta1, double tetta2, double m0, double h, double P, double Wisl, double hag, bool isNeedtoPrint)
{
    vector<double> t1;
    vector<double> t2;
    double eps = 0.1;
    double lymbda = 50;
    control Yat;
    t1.push_back(Mira.t1);
    t2.push_back(Mira.t2);
    int i = 0;
    Yat = kosnet(t1.back(), t2.back(), dtetta1, tetta2, m0, h, P, Wisl, hag);
    lol TT2 = rungekutta(t1.back(), t2.back(), dtetta1, tetta2, m0, h, P, Wisl, hag);
    do
    {
        double ddmdt1 = dmdt1(t1[i], t2[i], Yat.dtetta1, Yat.tetta2, m0, TT2.m, h, P, Wisl, hag, eps);
        double ddmdt2 = dmdt2(t1[i], t2[i], Yat.dtetta1, Yat.tetta2, m0, TT2.m, h, P, Wisl, hag, eps);
        cout << setprecision(7) << '\t' << m0 - TT2.m << '\t' << t1.back() << '\t' << t2.back() << '\t' << ddmdt1 << '\t' << ddmdt2 << '\n';
        t1.push_back(t1[i] - lymbda * ddmdt1);
        t2.push_back(t2[i] - lymbda * ddmdt2);
        if (t2.back() < t1.back())
        {
            cout << "Something went wrong" << '\n';
            system("pause");
        }
        Yat = kosnet(t1.back(), t2.back(), Yat.dtetta1, Yat.tetta2, m0, h, P, Wisl, hag);
        lol TT3 = rungekutta(t1.back(), t2.back(), Yat.dtetta1, Yat.tetta2, m0, h, P, Wisl, hag);
        double kubas = (abs(TT3.m - TT2.m));
        if ((kubas < 1e-3) && (sqrt(ddmdt1 * ddmdt1 + ddmdt2 * ddmdt2) < 1e-5))
        {
            cout << setprecision(7) << '\t' << m0 - TT2.m << '\t' << t1.back() << '\t' << t2.back() << '\t' << ddmdt1 << '\t' << ddmdt2 << '\n';
            TT3 = rungekutta(t1.back(), t2.back(), Yat.dtetta1, Yat.tetta2, m0, h, P, Wisl, hag, isNeedtoPrint);
            cout << setprecision(15) << TT3.m << '\n';
            return Yat;
        }
        TT2 = TT3;
        i++;
    } while (true);
}
double proizvm(double t1, double t2, double dtetta1, double tetta2, double m0, double h, double P, double Wisl, double hag)
{
    double dm = 0.1;
    control Yat, Mira;
    auto a = [&](control Mira, double m0)
    {
        lol TT = rungekutta(Mira.t1, Mira.t2, Mira.dtetta1, Mira.tetta2, m0, h, P, Wisl, hag);
        return TT.m;
    };
    Yat = kosnet(t1, t2, dtetta1, tetta2, m0 + dm, h, P, Wisl, hag);
    Mira = vremya(Yat, Yat.dtetta1, Yat.tetta2, m0 + dm, h, P, Wisl, hag);
    double mk1 = a(Mira, m0 + dm);
    Yat = kosnet(Mira.t1, Mira.t2, Mira.dtetta1, Mira.tetta2, m0 - dm, h, P, Wisl, hag);
    Mira = vremya(Yat, Yat.dtetta1, Yat.tetta2, m0 - dm, h, P, Wisl, hag);
    double mk2 = a(Mira, m0 - dm);
    double mp1 = mk1 - (Mcons + dm);
    double mp2 = mk2 - (Mcons - dm);
    return (mp1 - mp2) / (2 * dm);
}
double proizvWisl(double t1, double t2, double dtetta1, double tetta2, double m0, double h, double P, double Wisl, double hag)
{
    double dWisl = 0.1;
    control Yat, Mira;
    auto a = [&](control Mira, double Wisl)
    {
        lol TT = rungekutta(Mira.t1, Mira.t2, Mira.dtetta1, Mira.tetta2, m0, h, P, Wisl, hag);
        return TT.m;
    };
    Yat = kosnet(t1, t2, dtetta1, tetta2, m0, h, P, Wisl + dWisl, hag);
    Mira = vremya(Yat, Yat.dtetta1, Yat.tetta2, m0, h, P, Wisl + dWisl, hag);
    double mk1 = a(Mira, Wisl + dWisl);
    Yat = kosnet(Mira.t1, Mira.t2, Mira.dtetta1, Mira.tetta2, m0, h, P, Wisl - dWisl, hag);
    Mira = vremya(Yat, Yat.dtetta1, Yat.tetta2, m0, h, P, Wisl - dWisl, hag);
    double mk2 = a(Mira, Wisl - dWisl);
    double mp1 = mk1 - (Mcons);
    double mp2 = mk2 - (Mcons);
    return (mp1 - mp2) / (2 * dWisl);
}
#define _USE_MATH_DEFINES
#include <iostream>
#include "kepler.h"
#include <cmath>
#include <iomanip>
#include <fstream>
using namespace std;
// double mu = 398600.4418;// Земля
double mu = 4.903 * 1e3; // Луна
// double Rz = 6378.137;// Земля
double Rz = 1738;
double b0 = 398600.44;             // коэф 0 гармоники Земли
double b2 = 0.17552 * pow(10, 11); // коэф 2 гармоники Земли
decard keplertodec(double r, double omg, double u, double i, double TA, double e)
{
    TA = TA * M_PI / 180;   // истинная аномалия
    omg = omg * M_PI / 180; // долгота восходящего узла
    u = u * M_PI / 180;     // аргумент широты
    i = i * M_PI / 180;     // наклонение орбиты
    double a = r;
    double p = a * (1 - pow(e, 2));
    double Vr = pow(mu / p, 0.5) * e * sin(TA);
    double Vn = pow(mu * p, 0.5) / r;
    decard dekard;
    dekard.x = 0;
    dekard.y = 0;
    dekard.z = 0;
    dekard.Vx = 0;
    dekard.Vx = 0;
    dekard.Vz = 0;
    dekard.x = r * (cos(omg) * cos(u) - sin(omg) * sin(u) * cos(i));
    dekard.y = r * (sin(omg) * cos(u) + cos(omg) * sin(u) * cos(i));
    dekard.z = r * sin(u) * sin(i);
    dekard.Vx = Vr * (cos(omg) * cos(u) - sin(omg) * sin(u) * cos(i)) - Vn * (cos(omg) * sin(u) + sin(omg) * cos(u) * cos(i));
    dekard.Vy = Vr * (sin(omg) * cos(u) + cos(omg) * sin(u) * cos(i)) - Vn * (sin(omg) * sin(u) - cos(omg) * cos(u) * cos(i));
    dekard.Vz = Vr * sin(u) * sin(i) + Vn * cos(u) * sin(i);
    return dekard;
}
kepler dectokepler(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0)
{
    double r0 = sqrt(pow(x0, 2) + pow(y0, 2) + pow(z0, 2));
    double V0 = sqrt(pow(Vx0, 2) + pow(Vy0, 2) + pow(Vz0, 2));
    double k0 = r0 * pow(V0, 2) / mu;
    double tetta0 = asin((x0 * Vx0 + y0 * Vy0 + z0 * Vz0) / (r0 * V0)); // угол наклонения вектора скорости
    double sintetta0 = (x0 * Vx0 + y0 * Vy0 + z0 * Vz0) / (r0 * V0);
    double costetta0 = (sqrt(pow(r0, 2) * pow(V0, 2) - pow(x0 * Vx0 + y0 * Vy0 + z0 * Vz0, 2))) / (r0 * V0);
    if (tetta0 < -M_PI / 2 && tetta0 > M_PI / 2)
    {
        cout << "error";
    }
    kepler data;
    data.a = 0;   // большая полуось
    data.e = 0;   // эксцентриситет
    data.omg = 0; // долгота восходящего узла
    data.tau = 0; // время прохождения через перигей
    data.w = 0;   // аргумент перицентра
    data.TA = 0;  // истинная аномалия
    double u0 = 0;
    double TA = 0; // истинная аномалия
    data.a = r0 / (2 - k0);
    data.e = sqrt(pow(1 - k0, 2) + k0 * (2 - k0) * pow(sintetta0, 2));
    double C1 = y0 * Vz0 - z0 * Vy0;
    double C2 = z0 * Vx0 - x0 * Vz0;
    double C3 = x0 * Vy0 - y0 * Vx0;
    double C = sqrt(pow(C1, 2) + pow(C2, 2) + pow(C3, 2));
    data.i = acos(C3 / C);
    // дальнейшие условия объяснять нет смысла
    if (data.i < 0 && data.i > M_PI)
    {
        cout << "error";
    }
    if (abs(data.i - 0) < 1e-6)
    {
        data.omg = 0;
        data.w = 0;
        if (x0 < 0 || y0 < 0 || z0 < 0)
        {
            u0 = acos((x0 * cos(data.omg) + y0 * sin(data.omg)) / r0);
        }
        else
        {
            u0 = acos((x0 * cos(data.omg + M_PI) + y0 * sin(data.omg + M_PI)) / r0);
        }
        data.TA = u0;
    }
    else
    {
        if (x0 < 0 || y0 < 0 || z0 < 0 || abs(x0 - 0) < 1e-6 || abs(y0 - 0) < 1e-6 || abs(z0 - 0) < 1e-6)
        {
            data.omg = atan(-C1 / C2);
        }
        else
        {
            data.omg = atan(-C1 / C2) + M_PI;
        }
        if (data.e < 1e-6)
        {
            if (x0 < 0 || y0 < 0 || z0 < 0 || abs(x0 - 0) < 1e-6 || abs(y0 - 0) < 1e-6 || abs(z0 - 0) < 1e-6)
            {
                u0 = acos((x0 * cos(data.omg) + y0 * sin(data.omg)) / r0);
            }
            else
            {
                u0 = acos((x0 * cos(data.omg + M_PI) + y0 * sin(data.omg + M_PI)) / r0);
            }
            data.TA = u0;
            data.w = u0 - data.TA;
        }
        else
        {
            double sinTA = (k0 * sintetta0 * costetta0) / data.e;
            double cosTA = ((k0 * pow(costetta0, 2)) - 1) / data.e;
            double k8 = pow(sinTA, 2) + pow(cosTA, 2);
            data.TA = asin((k0 * sintetta0 * costetta0) / data.e);
            double TA1 = acos((k0 * pow(costetta0, 2) - 1) / data.e);
            double u0 = 0;
            if (x0 < 0 || y0 < 0 || z0 < 0 || abs(x0 - 0) < 1e-6 || abs(y0 - 0) < 1e-6 || abs(z0 - 0) < 1e-6)
            {
                double op = ((x0 * cos(data.omg) + y0 * sin(data.omg)) / r0);
                if (abs(op - 1) < 1e-6)
                {
                    op = 1;
                    u0 = acos(op);
                }
                else
                {
                    u0 = acos((x0 * cos(data.omg) + y0 * sin(data.omg)) / r0);
                }
            }
            else
            {
                u0 = acos((x0 * cos(data.omg + M_PI) + y0 * sin(data.omg + M_PI)) / r0);
            }
            data.w = u0 - data.TA;
        }
    }
    data.i = data.i * 180 / M_PI;
    if (data.omg * 180 / M_PI < 0)
    {
        data.omg = 360 + data.omg * 180 / M_PI;
    }
    else
    {
        data.omg *= 180 / M_PI;
    }
    data.w *= 180 / M_PI;
    if (data.TA * 180 / M_PI < 0)
    {
        data.TA = 360 + data.TA * 180 / M_PI;
    }
    else
    {
        data.TA *= 180 / M_PI;
    }
    return data;
}
pop Hpa(double a, double e)
{
    pop data;
    data.Hp = a * (1 - e);
    data.Ha = a * (1 + e);
    return data;
}
#ifndef kepler_h
#define kepler_h

class kepler
{
public:
    double a;
    double e;
    double w;
    double omg;
    double tau;
    double i;
    double TA;
};
class decard
{
public:
    double x;
    double y;
    double z;
    double Vx;
    double Vy;
    double Vz;
};
class Speed
{
public:
    double V1;
    double V2;
    double V;
};
class povort1
{
public:
    double V;
    double TA;
};
class BiSpeed
{
public:
    double V1;
    double V2;
    double V3;
};
class bob
{
public:
    double a;
    double e;
};
class pop
{
public:
    double Hp;
    double Ha;
};
decard keplertodec(double r, double omg, double u, double i, double TA, double e);
kepler dectokepler(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0);
double C(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double Hp1, double Ha1, double Hp2, double Ha2, double I, double m0, double i);
double povorot(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double Rp, double Ra, double i);
bob ae(double Hp2, double Ha2);
pop Hpa(double a, double e);
Speed Goman(double Hp1, double Ha1, double Hp2, double Ha2);
double gom(double Hp1, double Ha1, double Hp2, double Ha2, double Hp3, double Ha3, double i, double omg, double TA);
double V24(double TA, double e, double a, double i);
double dxdt(double Vx);
double dydt(double Vy);
double dzdt(double Vz);
double dVxdt(double x, double y, double z);
double dVydt(double x, double y, double z);
double dVzdt(double x, double y, double z);
void rungekutta(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double t0, double tk);
#endif // !
#pragma once

class lol
{
public:
    double r;
    double V;
    double stetta;
    double m;
};
class control
{
public:
    double t;
    double t1;
    double t2;
    double dtetta1;
    double tetta2;
};
class PA
{
public:
    double hagi;
    double dt;
};
double dVxdt(double t, double x, double y, double t1, double t2, double dtetta1, double tetta2, double m0, double P, double betta, bool j);
double dVydt(double t, double x, double y, double t1, double t2, double dtetta1, double tetta2, double m0, double P, double betta, bool j);
double fP(double, double, double, double P, bool);
double dxdt(double Vx);
double dydt(double Vy);
double dmdt(double m0, double betta, double t, double t1, double t2, bool j);
double gx(double x, double y);
double gy(double x, double y);
double tangag(double t, double t1, double t2, double dtetta1, double tetta2, bool);
double uglnaklonakstarkgor(double Vy, double Vx);
double Risl(double h);
double Visl(double h);
double fi(double Tetta, double uglnaklona);
double scalyartimes(double x, double y, double Vx, double Vy);
double time(double t, double hag, double t1, double t2, double i, double V, double Vorb);
lol rungekutta(double t1, double t2, double dtetta1, double tetta2, double m0, double h, double P, double Wisl, double hag, bool isNeedtoPrint = false);
control kosnet(double t1, double t2, double dtetta1, double tetta2, double m0, double h, double P, double Wisl, double hag);
control vremya(control Mira, double dtetta1, double tetta2, double m0, double h, double P, double Wisl, double hag, bool isNeedtoPrint = false);
double dmdt1(double t1, double t2, double dtetta1, double tetta2, double m0, double m, double h, double P, double Wisl, double hag, double eps);
double dmdt2(double t1, double t2, double dtetta1, double tetta2, double m0, double m, double h, double P, double Wisl, double hag, double eps);
bool flagept(double t, double t1, double t2);
double proizvm(double t1, double t2, double dtetta1, double tetta2, double m0, double h, double P, double Wisl, double hag);
double proizvWisl(double t1, double t2, double dtetta1, double tetta2, double m0, double h, double P, double Wisl, double hag);
double dVypravdt(double t, double t1, double t2, double m, double P, double alf);
double dVgdt(double x, double y, double Vx, double Vy);
