#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from math import sqrt

X = [32, 144, 63]
Y = [66, 182, 99]

# funkcja z dokumentają 
def odl3D(A, B) :
    """
    Opis Funkcji (1-2 zdania)
    
    Parameters:
    --------------
         A : Typ zmiennej - Opis [jednostka]
         B : Typ zmiennej - Opis [jednostka]
    
    Returns:
    --------------
        odl : Typ zmiennej - Opis [jednostka]

    """
    od = sqrt( (A[0] - B[0])**2 + (A[1] - B[1])**2 + (A[2] - B[2])**2 )
    return(od)

# wywołanie i przypisanie wyniku jej dziłania do zmiennej "odleglosc"
odleglosc = odl3D(X, Y)
print('wynik dziłania funkcji: ', odleglosc)
# DOSTĘP DO DOKUMENTACJI 
print(odl3D.__doc__) # 1 sposób dostępu do dokumentacji
print(help(odl3D))   # 2 sposób dostępu do dokumentacji



# funkcja z dokumentają  i adnotacja
def odl3D(A: 'cm', B: 'cm') ->  'cm' :
    """
    Funkcja do obliczenia długości na podstawie współrzędnych 3D.
    
    Parameters:
    --------------
         A : list of floats - współrzędne pkt A = [x, y, z] [cm]
         B : list of floats - współrzędne pkt B = [x, y, z] [cm]
    
    Returns:
    --------------
        odl - float - odległość AB [cm]

    """
    od = sqrt( (A[0] - B[0])**2 + (A[1] - B[1])**2 + (A[2] - B[2])**2 )
    return(od)


def decimalDeg2dms(decimalDeg):
    '''
    Metoda przelicza wartość kąta podaną w dziesiętnych stopnia 
    do wartości kata wyrażonej w stopniach, minutach i sekundach
    INPUT:
        decimalDeg [float]: decimal degree (np, 185.209060574544st)
    OUTPUT: 
        (d, m, s) [tuple]: (stopnie, minuty, sekundy) (np.(185st, 12min, 32.62sec))
    '''
    d = int(decimalDeg)
    m = int((decimalDeg - d) * 60)
    s = (decimalDeg - d - m / 60.) * 3600.0
    return (d, m, s) 


def dms2decimalDeg(st, minut, sec):
    '''
    Metoda przelicza wartość kąta podaną w w stopniach, minutach i sekundach 
    do wartości kata wyrażonej dziesiętnych stopnia 
    INPUT:
        st      : [int]   : wartość stopni
        minut   : [int]   : wartosc minut
        sec     : [float] : wartosc sekund 
    OUTPUT: 
        deg     : [float] : wartośc kąta w dziesiętnych stopnia
    '''
    deg = st + minut/60 + sec / 3600.
    return deg 


def azimuth_dist_xy(xA, yA, xB, yB):
    """
    Wyznaczenie azymutu AB i odległości skośniej pomiedzy punktami AB
    INPUT:
        xA : [float] : współrzędna x ounktu A
        yA : [float] : współrzędna y ounktu A
        xB : [float] : współrzędna x ounktu B
        yB : [float] : współrzędna y ounktu B
    OUTPUT:
        (Az_deg, dist_AB) - krotka dwuelementowa, gdzie:
            Az_deg : [float] : azymut AB w stopniach dziesiętnych
            dist_AB: [float] : odległość AB w jednostkach jak podano współrzędne.
    EXAMPLE:    
        INP: xA =-45.00; yA = 23.82; xB = 67.98; yB = 34.12 
        RUN: az, dist = azimuth_dist_xy(xA, yA, x_B, y_b)
        OUT: 5.209060574544288, 113.44853635018832
    """
    # wyznaczenie przyrostów współrzednych
    dX = xB - xA
    dY = yB - yA 
    # wyznaczenie azymutu:
    if dX > 0 and dY > 0:                   # I ćwiartka (0-90st)
        Az      = atan(dY/dX)               # [rad]
        Az_deg  = degrees(Az)               # [deg]
    elif dX < 0 and dY > 0:                 # II ćwiartka (90-180st)
        Az      = atan(dY/dX)+  pi          # [rad]
        Az_deg  = degrees(Az)               # [deg]
    elif dX < 0 and dY < 0:                 # III ćwiartka (180-270st)
        Az      =  atan(dY/dX) +  pi        # [rad]
        Az_deg  =  degrees(Az)              # [deg]
    elif dX > 0 and dY > 0:                 # IV ćwiartka (270-360st)
        Az      =  atan(dY/dX)  + 2 *  pi   # [rad]
        Az_deg  =  degrees(Az)              # [deg]
    elif dX == 0 and dY > 0:                # (90st)
        Az      =  pi /2                    # [rad]
        Az_deg  =  degrees(Az)              # [deg]
    elif dX < 0 and dY == 0:                # (180st)
        Az      =  pi                       # [rad]
        Az_deg  =  degrees(Az)              # [deg]
    elif dX == 0 and dY < 0:                # (270st)
        Az      =  pi +  pi /2              # [rad]
        Az_deg  =  degrees(Az)              # [deg]
    elif dX > 0 and dY == 0:                # (360st lub 0st)
        Az1     = 0                         # [rad]
        Az_deg1 =  degrees(Az)              # [deg]
        Az2     = 2*  pi                    # [rad]
        Az_deg2 =  degrees(Az)              # [deg]
    # wyznaczenie długości odcinka AB
    dist_AB =  sqrt(dX**2 +dY**2) # meter
    return Az_deg, dist_AB

def Np(phi, a, e2):
    """
    Promień krzywizny na pozycję uzytkownika
    Compute East-West Radius of curvature at current position
    INPUT:
        phi : [float] : szerokość geodezyjna (dziesiętne stopnia)
        a   : [float] : duża półoś elispoidy (promień równikowy) (metry)
        e2  : [float] : spłaszczenie elispoidy do kwadratu
    INPUT:
        N   : [float] : promien krzywizny
    """
    N = a/(1-e2*(sin(phi))**2)**(0.5)
    return N

def hirvonen(X, Y, Z, a = 6378137., e2 = 0.006694379990):
    """
    Algorytm Hirvonena - algorytm służący do transformacji współrzędnych ortokartezjańskich (prostokątnych) x, y, z 
    na współrzędne geodezyjne phi, lam, h. Jest to proces iteracyjny. 
    W wyniku 3-4-krotnego powtarzania procedury można przeliczyć współrzędne z dokładnoscią ok 1 cm.
 
    INPUT:
        X : [float] - współrzędna geocentryczna (ortokartezjański)
        Y : [float] - współrzędna geocentryczna (ortokartezjański)
        Z : [float] - współrzędna geocentryczna (ortokartezjański)
        a : [floar] - duża półoś elispoidy (promień równikowy) (metry)
        e2: [float] - spłaszczenie elipsoidy do kwadratu 
            
        inicjalizacji dla elipsoidy WGS84:
        a = 6378137
        e2= 0.0818191908426215**2 = 0.006694379990141318
    OUTPUT:
        phi :[float] : szerokość geodezyjna (stopnie dziesiętne)
        lab :[float] : długość geodezyjna (stopnie dziesiętne)
        hel :[float] : wysokość elipsoidalna (metry)
    EXAMPLE: 
        INP: X = 3731440.0; Y = 1240560.0; Z = 5005620.0 
        RUN: phi, lam, H = hirvonen(X, Y, Z) 
        OUT: 52.034738613586406, 18.389978007504855, 555.4033404150978
        
    CONTROL:
        jeśli sqrt(X**2 + Y**2 + Z**2)
    """
    r   = sqrt(X**2 + Y**2)           # promień
    phi = atan(Z / (r * (1 - e2)))    # pierwsze przyphilizenie
    j=0
    phi_next = phi
    while True:
        j+=1
        phi_prev = phi_next
        N    = Np(phi_prev, a, e2)
        hel  = (r/cos(phi_prev))- N
        phi_next = atan(Z/(r *(1 - (e2 * (N/(N + hel))))))
        phi_next = phi_next
        if abs(phi_next - phi_prev) < (0.0000001/206265):  # rho'' =206265
            break
    phi = phi_prev
    lam   = atan(Y/X)
    N   = Np(phi, a, e2)
    hel = (r/cos(phi))- N
    return degrees(phi), degrees(lam), hel



if __name__=='__main__':
    '''
     __name__, przechowuje nazwę modułu (wykonanego pliku pythona)
    Jesli wykonujemy (kompilujemy) ten plik to jego nazwa, nazwa modułu jest '__main__', więc:
    warunek "if __name__=='__main__':" >>> True
    '''
    wynik1 = decimalDeg2dms(185.20906057454428) # wynik: (185, 12, 32.618068359413385)
    print('wynik1 z main', wynik1)

    wynik2 = dms2decimalDeg(wynik1[0], wynik1[1], wynik1[2])
    print('wynik2 z main', wynik2)
    
    xyz = [3731440.0,  1240560.0, 5005620.0]
    phi, lam, hel = hirvonen(xyz[0], xyz[1], xyz[2], a = 6378137., e2 = 0.006694379990)

    # wywołanie funkcji i przypisanie wyniku jej dziłania do zmiennej odleglosc
    odleglosc = odl3D(X, Y)

    # DOSTĘP DO ADNOTACJI 
    adnotacje = odl3D.__annotations__ # wywołanie słownika z adnotacjami
    print('Adnotacje funkcji sa dostępne w słowniku: \n',  adnotacje)
    print('Odnosząc się do własciwego klucza uzyskamy adnotacje \n',  adnotacje['return'])

    # w ten sposób można np. podać wynik z jednostką
    print('Odległośc pomiędzy punktem X i Y wynosi:', round(odleglosc, 3),  adnotacje['return'])

