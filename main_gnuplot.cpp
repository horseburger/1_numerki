/*
Niniejszy program znajduje miejsca zerowe dla różnych funkcji metodą styczncych i bisekcji.

Wykonane przez :
	Jakub Zygmunt
	Stefan Żaryń
	FTIMS INFORMATYKA
*/

#include <iostream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits.h>
#include "gnuplot_i.hpp"

#define DOKLADNOSC 0
#define ITERACJE 1

using namespace std;

double a,b,eps; 		//przedzialy
int iter;				//iteracje
string funkcje [4]={
	"x*x-4*x",
	"sin(x)",
	"(2**x)-5",
	"sin(x)*sin(x)-4*sin(x)"
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

void menu(void); //wyswietla menu
void wykres(double a,double b, double zer0, double(*f)(double), string nazwa); //rysuje wykres w gnuplot dla wartosci podanych przez uzytkownika
void wykres(string nazwa); //rysuje wykres w gnuplot odpowiedniej funkcji aby uzytkownik mogl znalezc przedzial gdzie jest miejsce zerowe
void select(int & wybor); //switch dla wybranej funkcji przez uzytkownika
void metoda (double a,double b, double(*f)(double), double(*fp)(double), string nazwa); //uzytkownik wybiera czy chce podac wartosc epsilon czy ilosc iteracji
void podsumowanie(const double & a, const double & b); //krotkie porownanie obu metod

//bisekcja
double w_bisekcja (double a, double b, double( * f )( double ),int wybor, double eps = -1, unsigned int iter = UINT_MAX);

//styczne
double w_styczne (double a,double b, double( * f )( double ), double ( * fp )( double ),int wybor, double eps=-1, unsigned int iter=UINT_MAX);

int main()
{
	//
	cout << setprecision(8) << fixed;
	int wybor;
	char koniec;
	cout << "Wpisz 'k' aby zakonczyc program lub dowolny klawisz aby kontynuowac : ";
	while(cin.get(koniec) && koniec != 'k')
	{
	menu();
	cin >> wybor;
	getchar(); // Usunięcie znaku końca linii zostawionego przez cin
	select(wybor);
	cout << "\nWpisz 'k' aby zakonczyc program lub dowolny klawisz aby kontynuowac : ";
	}
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
			double bis_eps = w_bisekcja(a, b, f, 0, eps);
			double st_eps  = w_styczne(a, b, f, fp, 0,eps);
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
			double bis_it = w_bisekcja(a, b, f, 1, -1, iter);
			double st_it  = w_styczne(a, b, f, fp, 1,-1,iter);
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
			 << "4. Wielomian(trygonometryczna)" << endl;
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
			default:
				cout << "Bledna wartosc\n";
				break;
		}
}
double w_bisekcja (double a, double b, double( * f )( double ),int wybor, double eps, unsigned int iter)
{
	double wymagana_dokladnosc = -1;
	unsigned int maks_iteracji = UINT_MAX;
	if ( DOKLADNOSC == wybor )
	{
		wymagana_dokladnosc = eps;
	} else if ( ITERACJE == wybor ) {
		maks_iteracji = iter;
	}
	double fa, fb, f0;
	double zerowa;
		fa = f(a);
		fb = f(b);
	double aktualna_dokladnosc = fabs(a - b);
	int ilosc_iteracji = 0;

	 while ( ( ilosc_iteracji < maks_iteracji ) && ( aktualna_dokladnosc > wymagana_dokladnosc ) )
			{
			ilosc_iteracji++;
			zerowa = (a + b) / 2.0;
			f0 = f(zerowa);
			if(fabs(f0) <= wymagana_dokladnosc)
			{
			cout << "Wykonalem metoda bisekcji " << ilosc_iteracji << " iteracji\n";
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
			cout << "Wykonalem metoda bisekcji " << ilosc_iteracji << " iteracji\n";
			return zerowa;


}
double w_styczne (double a,double b, double( * f )( double ), double ( * fp )( double ),int wybor, double eps, unsigned int iter)
{
		double wymagana_dokladnosc = -1;
		unsigned int maks_iteracji = UINT_MAX;
		if ( DOKLADNOSC == wybor )
		{
			wymagana_dokladnosc = eps;
		} else if ( ITERACJE == wybor ) {
			maks_iteracji = iter;
		}
		double x0,x1,f0,f1;
		x0=a;
		x1 = x0 - 1;
		f0 = f(x0);
		unsigned int ilosc_iteracji = 0;

		 while ( ( ilosc_iteracji < maks_iteracji ) && ( fabs(f0) > wymagana_dokladnosc ) && (fabs(x1 - x0) > wymagana_dokladnosc ))
				{
						ilosc_iteracji++;
						f1 = fp(x0);
						if(fabs(f1) < wymagana_dokladnosc)
						{
							cout << "Zly punkt startowy dla metody stycznych";
							ilosc_iteracji = UINT_MAX;
						}
						x1 = x0;
						x0 = x0 - (f0 / f1);
						f0 = f(x0);
				}
				cout << "Wykonalem metoda stycznych " << ilosc_iteracji << " iteracji\n";
				return x0;

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