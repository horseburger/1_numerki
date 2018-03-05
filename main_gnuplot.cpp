#include <iostream>
#include <cmath>
#include <cstring>
#include "gnuplot_i.hpp"

using namespace std;

double a,b,eps; 		//przedzialy
int iter;				//iteracje	
string funkcje [9]={
	"x*x-4*x",
	"sin(x)",
	"(2**x)-5",
	"sin(x)*sin(x)-4*sin(x)",
	"((2**x)-5)*((2**x)-5)-4*((2**x)-5)",
	"sin(x*x-4*x)",
	"sin((2**x)-5)",
	"(2**(x*x-4*x))-5",
	"(2**sin(x))-5"
};						//funkcje

//funkcje + jej pochodne
double fun_nieliniowa(double x);
double p_fun_nieliniowa(double x);

double trygonmetryczna(double x);
double p_trygonmetryczna(double x);

double wykladnicza(double x);
double p_wykladnicza(double x);

double nieliniowa_t(double x);
double p_nieliniowa_t(double x);

double nieliniowa_w(double x);
double p_nieliniowa_w(double x);

double trygonometryczna_n(double x);
double p_trygonometryczna_n(double x);

double trygonometryczna_w(double x);
double p_trygonometryczna_w(double x);

double wykladnicza_n(double x);
double p_wykladnicza_n(double x);

double wykladnicza_t(double x);
double p_wykladnicza_t(double x);



void menu(void);
void sub_menu();
void wykres(double a,double b, double zer0, double(*f)(double), string nazwa);
void wykres(string nazwa);
void wykres(double a, double b, double zer0, double(*f)(double),double(*g)(double),string nazwa);
void select(int & wybor);
void metoda (double a,double b, double(*f)(double), double(*fp)(double), string nazwa);
void podsumowanie(const double & a, const double & b);

//bisekcja
double bisekcja (double a, double b, double eps, double( * f )( double ));
double bisekcja (double a, double b, int iter, double( * f )( double ));

//styczne
double styczne (double a, double eps, double( * f )( double ), double( * fp )( double ));
double styczne (double a, int iter, double( * f )( double ), double( * fp )( double ));

int main() 
{
	//
	
	int wybor;
	menu();
	cin >> wybor;
	getchar(); // Usunięcie znaku końca linii zostawionego przez cin
	select(wybor);
	
	//
return 0;	
}
void metoda (double a,double b, double(*f)(double), double(*fp)(double), string nazwa)
{
	if((f(a)*f(b)) <= 0)
	{
			cout << "Jaka metoda wyznaczania miejsca zerowego ?\n";
			cout << "1.Przez podanie wartosci epsilon\n";
			cout << "2.Przez ilosc iteracji\n";
			int metoda;
			cin >> metoda;

		if(metoda == 1)
		{
			cout << "Wartosc epsilon (zblizony do 0, np 0.0001) : ";
			cin >> eps;
			getchar();
			double bis_eps = bisekcja(a,b,double(eps),f);
			double st_eps  = styczne(a, double(eps), f, fp);
			if(eps)
			{
			wykres(a,b,bis_eps,f,nazwa+(" bisekcja"));
			wykres(a, b,st_eps , f, nazwa+(" styczne"));
			podsumowanie(bis_eps, st_eps);
			}
		}
		else if(metoda == 2)
		{
			cout << "Ilosc powtorzen : ";
			cin >> iter;
			getchar();
			double bis_it = bisekcja(a,b,iter,f);
			double st_it  = styczne(a, iter, f, fp);
			wykres(a, b, bis_it,f,nazwa+(" bisekcja"));
			wykres(a, b, st_it, f, nazwa+(" styczne"));
			podsumowanie(bis_it, st_it);
		}
		else
		{
		cout << "Zła wartosc!\n";
		}
	}
	
	else
	{
		cout << "Funkcja nie spelnia zalozen\n";
	}

}


void menu(void)
{
	cout << "Wybierz rodzaj funkcji: " << endl
			 << "1. Wielomian" << endl
			 << "2. Trygnometryczna" << endl
			 << "3. Wykladnicza" << endl
			 << "4. Wielomian(trygonometryczna)" << endl
			 << "5. Wielomian(wykladnicza)" << endl
			 << "6. Trygonometryczna(wielomian)" << endl
			 << "7. Trygonometryczna(wykladnicza)" << endl
			 << "8. Wykladnicza(wielomian)" << endl
			 << "9. Wykladnicza(trygonometryczna)" << endl;
}


void wykres(double a,double b,double zer0, double(*f)(double),string nazwa)
{
	Gnuplot wyk;
	wyk.set_grid();
	wyk.set_style("lines");
	wyk.set_xrange(-10,10);
	wyk.set_yrange(-10,10);
	wyk.set_title(nazwa);
	wyk.set_xlabel("Os x");
	wyk.set_ylabel("Os y");

	double tmp = a,roz = (abs(a-b)*5)+1;
	vector<double> x(roz);
	for(int i=0;i<roz;i++)
	{
		x[i]=tmp;
		tmp+=0.2;
	}

	vector<double> y(roz);
	for(int i=0;i<roz;i++)
	{
		y[i] = f(x[i]);
	}

	wyk.plot_xy(x,y,nazwa);

	vector<double> zer_x(1);
	vector<double> zer_y(1);
	zer_y[0] = 0;
	zer_x[0] = zer0;
	
	wyk.set_style("points");
	wyk.set_pointsize(5.0);
	wyk.plot_xy(zer_x,zer_y,"Miejsca zerowe");
	
	getchar();
}

void wykres(string nazwa)
{
	Gnuplot wyk;
	wyk.set_grid();
	wyk.set_style("lines");
	wyk.set_xrange(-10,10);
	wyk.set_yrange(-10,10);
	wyk.set_title(nazwa);
	wyk.set_xlabel("Os x");
	wyk.set_ylabel("Os y");
	wyk.plot_equation(nazwa);
	getchar();
}


void select(int & wybor)
{
	wykres(funkcje[wybor-1]);
	cout << "Podaj lewą strone przedziału : ";
	cin >> a;
	cout << "Podaj prawa strone przedziału : ";
	cin >> b;
	getchar();

	switch(wybor)
		{
			case 1:
			{
				metoda(a, b, fun_nieliniowa, p_fun_nieliniowa, funkcje[wybor-1]);
				break;
			}
			case 2:
			{
				metoda(a, b, trygonmetryczna, p_trygonmetryczna, funkcje[wybor-1]);
				break;
			}
			case 3:
			{
				metoda(a, b, wykladnicza, p_wykladnicza, funkcje[wybor-1]);
				break;
			}
			case 4:
			{
				metoda(a, b, nieliniowa_t, p_nieliniowa_t, funkcje[wybor-1]);
				break;
			}
			case 5:
			{
				metoda(a, b, nieliniowa_w, p_nieliniowa_w, funkcje[wybor-1]);
				break;
			}
			case 6:
			{
				metoda(a, b, trygonometryczna_n, p_trygonometryczna_n, funkcje[wybor-1]);
				break;
			}
			case 7:
			{
				metoda(a, b, trygonometryczna_w, p_trygonometryczna_w,funkcje[wybor-1]);
				break;
			}
			case 8:
			{
				metoda(a,b, wykladnicza_n, p_wykladnicza_n, funkcje[wybor-1]);
				break;
			}
			case 9:
			{
				metoda(a,b, wykladnicza_t, p_wykladnicza_t, funkcje[wybor-1]);
				break;
			}
			default:
				cout << "Bledna wartosc\n";
				break;
		}
}

double bisekcja (double a, double b, double eps, double( * f )( double ))
{
	int i=0;
	double fa, fb, f0;
	double zerowa;
	fa = f(a);
	fb = f(b);
		while(fabs(a - b) > eps)
			{
				i++;
				zerowa = (a + b) / 2.0; 
				f0 = f(zerowa);
				if(fabs(f0) < eps)
				{
				cout << "Wykonalem metoda bisekcji " << i << " iteracji\n";
				return zerowa;
				}
				if(fa * f0 < 0)
				b = zerowa;
				else
				{
					a = zerowa;
					fa = f0;
				}
			}
	
	cout << "Wykonalem metoda bisekcji " << i << " iteracji\n";
	return zerowa;
}	
double bisekcja (double a, double b, int iter, double( * f )( double ))
{
		double fa, fb, f0;
		double zerowa;
		fa = f(a);
		fb = f(b);
			while( (iter--) > 0)
				{
					zerowa = (a + b) / 2.0; 
					f0 = f(zerowa);
					if(fabs(f0) == 0)
					return zerowa;
					if(fa * f0 < 0)
					b = zerowa;
					else
					{
						a = zerowa;
						fa = f0;
					}
				}
			return zerowa;
}
double styczne (double a, double eps, double( * f )( double ), double ( * fp )( double ))
{
	int i =0;
	double x, f0, f1;
	x = a - 1; 
	f0 = f(a); 
		while ((fabs(x - a) > eps) && (fabs(f0) > eps))
		{
			f1 = fp(a);
			if(fabs(f1) < eps)
			{
				cout << "Zly punkt startowy\n";
				break;
			}
			x = a;
			a = a - (f0 / f1);
			f0 = f(a);
			i++;
		}
		cout << "Wykonalem metoda stycznych " << i << " iteracji\n";
	return x;
}
double styczne (double a, int iter, double( * f )( double ), double( * fp )( double ))
{
	double x, f0, f1;
	x = a - 1; 
	f0 = f(a); 
		while (iter--)
		{
			f1 = fp(a);
			if(fabs(f1) == 0)
			{
				cout << "Zly punkt startowy\n";
				break;
			}
			x = a;
			a = a - (f0 / f1);
			f0 = f(a);
		}
	return x;
}
void podsumowanie(const double & a, const double & b)
{
	cout << "Miejsce zerowe dla metody bisekcji : " << a << endl ;
	cout << "Miejsce zerowe dla metody stycznych : " << b << endl ;
	double x = fabs(a-b);
	cout << "Roznica wynikow metod bisekcji i stycznych to : " << x << endl ;	
}
double fun_nieliniowa(double x)
{
	return (x*x) - (4*x);
}
double p_fun_nieliniowa(double x)
{
	return 2*x - 4;
}
double trygonmetryczna(double x)
{
	return sin(x);
}
double p_trygonmetryczna(double x)
{
	return cos(x);
}
double wykladnicza(double x)
{
	return (pow(2,x)-5);
}
double p_wykladnicza(double x)
{
	return pow(2,x)*log(2);
}
double nieliniowa_t(double x)
{
	return sin(x)*sin(x)-4*sin(x);
}
double p_nieliniowa_t(double x)
{
	return 2*(sin(x)-4)*cos(x);
}
double nieliniowa_w(double x)
{
	return ((pow(2,x)-5)*(pow(2,x)-5)-4*(pow(2,x)-5));
}
double p_nieliniowa_w(double x)
{
	return (pow(2,x+1)*(pow(2,x)-2)*log(2));
}
double trygonometryczna_n(double x)
{
	return sin(x*x-4*x);
}
double p_trygonometryczna_n(double x)
{
	return (2*(x-2)*cos(x*(x-4)));
}
double trygonometryczna_w(double x)
{
	return sin(pow(2,x)-5);
}
double p_trygonometryczna_w(double x)
{
	return pow(2,x)*log(2)*cos(pow(2,x)-5);
}
double wykladnicza_n(double x)
{
	return pow(2,x*x-4*x);
}
double p_wykladnicza_n(double x)
{
	return pow(2,x*x-4*x+1)*log(2)*(x-2);
}
double wykladnicza_t(double x)
{
	return pow(2,sin(x));
}
double p_wykladnicza_t(double x)
{
	return log(2)*pow(2,sin(x))*cos(x);
}

	