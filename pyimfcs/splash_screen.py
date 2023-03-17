#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 10:50:50 2023

@author: aurelienb
"""

import os
import sys

from PyQt5.QtWidgets import QGridLayout, QLabel, QDialog, QApplication, QSplashScreen

from PyQt5.QtCore import Qt
from PyQt5.QtGui import  QPixmap

art0 = """
    ZZZzz
      |\      _,,,---,,_
      /,`.-'`'    -.  ;-;;,_
     |,4-  ) )-,_. ,\ (  `'-'
    '---''(_/--'  `-'\_)  Felix Lee 
"""
art1=r"""
           __..--''``---....___   _..._    __
 /// //_.-'    .-/";  `        ``<._  ``.''_ `. / // /
///_.-' _..--.'_    \                    `( ) ) // //
/ (_..-' // (< _     ;_..__               ; `' / ///
 / // // //  `-._,_)' // / ``--...____..-' /// / //
"""

art2=r"""
           .               ,.
          T."-._..---.._,-"/|
          l|"-.  _.v._   (" |
          [l /.'_ \; _~"-.`-t
          Y " _(o} _{o)._ ^.|
          j  T  ,--.  T  ]
          \  l ( /-^-\ ) !  !
           \. \.  "~"  ./  /c-..,__
             ^r- .._ .- .-"  `- .  ~"--.
              > \.                      \
              ]   ^.                     \
              3  .  ">            .       Y  -Row
 ,.__.--._   _j   \ ~   .         ;       |
(    ~"-._~"^._\   ^.    ^._      I     . l
 "-._ ___ ~"-,_7    .Z-._   7"   Y      ;  \        _
    /"   "~-(r r  _/_--._~-/    /      /,.--^-._   / Y
    "-._    '"~~~>-._~]>--^---./____,.^~        ^.^  !
        ~--._    '   Y---.                        \./
             ~~--._  l_   )                        \
                   ~-._~~~---._,____..---           \
                       ~----"~       \
                                      \
"""
art_credits = "<a href=https://www.asciiart.eu/animals/cats> https://www.asciiart.eu/animals/cats  </a>"
class LoadingWindow(QSplashScreen):
    """Opens a splash screen to let the user know that the program is running"""
    def __init__(self):
        super().__init__()
        
        self.textLabel = QLabel("Processing, please wait")
        self.logo = QLabel(art1)
        self.logo.setTextFormat(0)
        self.art_credits =QLabel(art_credits)
        # set the layout
        layout = QGridLayout()        
        self.setLayout(layout)
        layout.addWidget(self.textLabel,0,0)
        layout.addWidget(self.logo,1,0)
        layout.addWidget(self.art_credits,2,0)

if __name__=="__main__":
    app = QApplication([])
    win = LoadingWindow()
    win.show()
    app.exec_()