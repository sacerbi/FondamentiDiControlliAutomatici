# Documentazione - Cartella Capitolo 3 (C++)

## Prerequisiti

Questo progetto utilizza **Eigen**, una libreria C++ per l'algebra lineare. Per ulteriori informazioni e per scaricare la libreria, consultare il [sito ufficiale di Eigen](https://eigen.tuxfamily.org/).

Per la stampa dei grafici si usa gnuplot.

## Installazione di Eigen

Assicurarsi di avere la cartella `Eigen` disponibile nel vostro sistema. È possibile scaricarla dal repository ufficiale e posizionarla nel percorso desiderato.

## Compilazione

Per compilare i file C++ in questa cartella, utilizzare il seguente comando generalizzato:
```bash
g++ ./<filename>.cpp -o <output_name> -I./Eigen/ -std=c++20
```

Dove:
- `<filename>` è il nome del file sorgente C++ (es. `Esempio3_3.cpp`)
- `<output_name>` è il nome dell'eseguibile desiderato (es. `Esempio3_3`)