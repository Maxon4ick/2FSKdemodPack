# 2FSKdemodPack
В данном задании было необходимо реализовать демодулятор 2FSK сигнала при пакетной передачи полезного сигнала. Старт посылка является нулем в начале сообщения,
стоп посылкой являются полторы единичной посылки. Заданный сигнал обладает достаточно большим отношением сигнал-шум, однако канал передачи имеет замирания, механизм
противодействия которым также реализован в программе. 
Для корректной работы тестов, необходимо добавить тестовые файлы из папки TestFiles в файлы собранного проекта.
При реализации демодулятора была использована открытая библиотека, написанная на языке С, fftw3lib
