# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 00:12:11 2021

В модуле реализованы функции, описывающие радиальную структуру диска в 
соответствии с моделью Дудорова и Хайбрахманова.

См. статью Dudorov, Khaibrakhmanov, 2014, Ap&SS, 352, 103.

@author: Sergey Khaibrakhmanov
"""

import numpy as np
import sys
sys.path.append("./")
# подключение модуля с физическими константами
from const import G, k, m_p, M_sun, au, mu, Rg

class Disk:
    def __init__(self, m, mdot, alpha, l_star):
        """
        Конструктор класса

        Parameters
        ----------
        m : число
            масса звезды, в массах Солнца.
        mdot : число
            масса звезды, в массах Солнца.
        alpha : число
            масса звезды, в массах Солнца.
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
        self.C_Teff_v = 150.0 * (self.m_dot * self.m)**0.25
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