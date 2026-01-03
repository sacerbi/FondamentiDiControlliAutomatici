# Istruzioni per la compilazione

Per compilare il programma è sufficiente il comando:
<code>g++ EsercizioX_Y.cpp -o EsercizioX_Y -I/usr/include/python3.12 -I/usr/lib/python3/dist-packages/numpy/core/include -lpython3.12 -std=c++17 </code>

In caso di errori, verificare i requisiti:
* <code>sudo apt-get install python3-numpy</code>
* <code>sudo apt-get install python3-matplotlib</code>

Il path di inclusione può variare a seconda del proprio set-up. Nella compilazione bisogna includere il path di Python, necessario a Matplotlib-cpp e il path di numpy.

Per individuare il path di numpy, ad esempio, si può eseguire il seguente comando python:
<code>python3 -c "import numpy; print(numpy.get_include())"</code>