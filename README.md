# Proiectare Filtre IIR și FIR prin Metode de Transformare

Proiect realizat pentru disciplina **Prelucrarea Semnalelor (PS)**, în cadrul Facultății de Automatică și Calculatoare, Universitatea Națională de Știință și Tehnologie POLITEHNICA București.

## Obiectivele Proiectului
* **Proiectarea filtrelor IIR** utilizând prototipuri analogice clasice (Butterworth, Cebîșev Tip I/II, Cauer/Eliptic) prin tehnica transformării biliniare (Tustin).
* **Studiu comparativ avansat** între structurile IIR și filtrele FIR (proiectate prin metoda ferestrei Hamming și metoda celor mai mici pătrate - CMMP) din punct de vedere al ordinului minim și al liniarității fazei.
* **Analiza dependenței** funcției de transfer discrete în raport cu variația toleranțelor și a perioadei de eșantionare $T_s$.
* **Optimizare multi-criterială** printr-un algoritm iterativ de căutare a ordinului minim, evaluat cu o funcție de cost complexă.

## Structura Codului & Implementare

Proiectul este structurat pe faze de dezvoltare, conținând scripturi directoare (Main-uri) și rutine matematice dedicate:

### Scripturi Principale de Simulare (Main-uri `/src`)
* `Faza1.m` — Gestionează **Faza 1** a proiectului. Calculează parametrii analogici predistorsionați pentru filtrul trece-jos etalon Butterworth, apelează funcția de proiectare, verifică grafic gabaritul și realizează prima comparație de ordin și eroare de spectru/fază cu metodele FIR (Fereastră și CMMP).
* `Faza2.m` — Coordonează **Faza 2** (Filtre Cauer și Cebîșev). Rulează simulările pentru specificațiile modificate (toleranța în banda de stopare $\Delta_s$ dublată) și generează automat matrici grafice de comparație (2x3 spectre de amplitudine în dB și 2x3 răspunsuri în fază) salvate în directorul de figuri.
* `Faza3.m` — Reprezintă **Concursul de Proiectare**. Definește noile specificații strânse de gabarit ($\Delta\omega \approx 0.095$ rad), apelează iterativ toate rutinele de căutare a ordinului minim, calculează costul fiecărui filtru, sortează structurile de la cea mai eficientă la cea mai slabă și complotează spectrul comparativ global.
* `Faza4.m` — Scriptul de test pentru **Suplimentul Teoretic**. Simulează comportamentul unui filtru Butterworth modificat pentru a asigura un câștig staționar neunitar în curent continuu ($|H(0)| = 1 + \Delta_p$).

### Algoritmi Iterativi de Căutare și Validare (`/src`)
Pentru a garanta respectarea strictă a toleranțelor impuse prin testarea dinamică a răspunsului în frecvență (`freqz`), s-au implementat următoarele funcții:
* `But_FTI.m` & `But_FTI_v2.m` — Soluția analitică Butterworth folosind transformarea Tustin standard, respectiv varianta modificată (fără factorul 2).
* `cauta_ordin_min_Cheby1.m` / `cauta_ordin_min_Cheby2.m` — Sinteza filtrelor Cebîșev prin verificarea automată a benzii complementare celei asigurate nativ de funcțiile MATLAB.
* `cauta_ordin_min_Eliptic.m` — Determinarea ordinului minim pentru filtrele Cauer (`ellip`), verificând dacă atenuarea $R_s$ este atinsă exact la frecvența limită a benzii de stopare.
* `cauta_ordin_min_Fereastra.m` / `cauta_ordin_min_CMMP.m` — Algoritmi iterativi pentru structuri FIR (ferestre Hamming `fir1` și optimizări Least Squares `firls`).
* `But_FTI_Faza4.m` — Modificarea formulelor de calcul pentru polii analogici și a pulsatației de tăiere $\Omega_c$ pentru integrarea amplificării de curent continuu direct în topologia filtrului.

### Funcția de Evaluare a Performanței
* `calc_criteriu_performanta.m` — Funcție centrală care calculează scorul de penalizare (Costul) în cadrul concursului, ponderând complexitatea hardware (60%), eroarea RMS de magnitudine în banda de trecere (30%) și eroarea RMS de fază față de o caracteristică ideal liniară (10%):
    $$C = 0.6 \cdot M + 0.3 \cdot (E_{mag}) + 0.1 \cdot (E_{phs})$$

## Documentație Completă
Demonstrațiile matematice complete pentru sistemul de ecuații nelineare de la Faza 4, tabelele centralizatoare de date, interpretările fenomenelor fizice (cum ar fi instabilitatea filtrelor IIR de ordin mare) și concluziile inginerești se găsesc în raportul atașat:
 **[Deschide Documentație PDF](Proiectarea_filtrelor_IIR_prin_metode_de_transformare-Martin_Razvan-Stefan.pdf)**

## Tehnologii Utilizate
* MATLAB R2024b (Signal Processing Toolbox)