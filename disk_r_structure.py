# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 00:12:11 2021

В модуле реализованы функции, описывающие радиальную структуру диска в 
соответствии с моделью Дудорова и Хайбрахманова.

См. статью Dudorov, Khaibrakhmanov, 2014, Ap&SS, 352, 103.

@author: Sergey Khaibrakhmanov
"""

import numpy as np
from scipy.interpolate import interp1d

class Disk:
    def __init__(self, m, mdot, alpha, l_star):
        """
        Конструктор класса

        Parameters
        ----------
        m : число
            масса звезды, в массах Солнца.
        mdot : число
            темп аккреции, в единицах 10^(-8)*Msun/год.
        alpha : число
            параметр турбулетности, в единицах 0.01.
        l_star : число
            светимость звезды, в светимостях Солнца.
            
        Returns
        -------
        None.

        """
        self.m = m         
        self.mdot = mdot   
        self.alpha = alpha 
        self.l_star = l_star
        
    def calc_const(self):
        """
        Вычисление вспомогательных констант в аналитическом решении уравнений модели диска
        Returns
        -------
        None.

        """
        
        self.C_T = 240 * self.m**(0.375) * self.alpha**(-0.25) * self.mdot**(0.5)
        self.C_Teff_v = 150.0 * (self.mdot * self.m)**0.25
        self.C_Teff_irr = 280.0 * self.l_star**0.25
        self.C_H = 0.03 * self.alpha**(-0.125) * self.mdot**(0.25) * self.m**(-0.3125)
    
    def set_model_params(self, m_in, mdot_in, alpha_in):
        """
        Установка параметров модели

        Parameters
        ----------
        m_in : число
            масса звезды, в массах Солнца.
        mdot_in : число
            темп аккреции, в единицах 10^(-8)*Msun/год.
        alpha_in : число
            турбулентный параметр Шакуры и Сюняева.

        Returns
        -------
        None.

        """
        self.m = m_in
        self.m_dot = mdot_in
        self.alpha = alpha_in

    def import_data(self, file_name):
        """
        Чтение файла данных с результатами расчетов структуры диска

        Parameters
        ----------
        file_name : строка
            имя файла данных.

        Returns
        -------
        None.

        """
        data = np.loadtxt(file_name, skiprows=1)
        self.r = data[:, 0]
        Teff = data[:, 6]
        
        method = 'cubic'
        self.T_i = interp1d(self.r, Teff, method)
        
    # ----- Радиальные профили основных величин -----
    
    def T(self, r_au):
        """
        Зависимость температуры газа от радиального расстояния
    
        Parameters
        ----------
        r_au : число
            радиальное расстояние от звезды, в а.е.
    
        Returns
        -------
        Температура газа на заданном расстоянии r_au, К.
    
        """
        return self.C_T * r_au**(-1.125)
    
    def Teff_v(self, r_au):
        """
        Зависимость эффективной температуры газа от радиального расстояния для случая вязкого нагрева
    
        Parameters
        ----------
        r_au : число
            радиальное расстояние от звезды, в а.е.
    
        Returns
        -------
        Температура газа на заданном расстоянии r_au, К.
    
        """
        return self.C_Teff_v * r_au**(-0.75)
    
    def Teff_irr(self, r_au):
        """
        Зависимость эффективной температуры газа от радиального расстояния для случая нагрева диска излучением звезды
    
        Parameters
        ----------
        r_au : число
            радиальное расстояние от звезды, в а.е.
    
        Returns
        -------
        Температура газа на заданном расстоянии r_au, К.
    
        """
        return self.C_Teff_irr * r_au**(-0.5)

    def Tnum(self, r_au):
        """
        Зависимость эффективной температуры газа от радиального расстояния из численного расчета структуры диска
    
        Parameters
        ----------
        r_au : число
            радиальное расстояние от звезды, в а.е.
    
        Returns
        -------
        Температура газа на заданном расстоянии r_au, К.
    
        """
        return self.T_i(r_au)
    
    def Get_r_in(self):
        return self.r[0]
    
    def Get_r_out(self):
        return self.r[(self.r).size - 1]