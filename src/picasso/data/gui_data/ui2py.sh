#!/bin/bash
#
# SCRIPT TO INSERT THIS AFTER UI to PY CONVERSION:
# import pkg_resources
# png_of_alpha = pkg_resources.resource_filename('picasso.src.picasso.data', 'gui_data/of_alpha.png')
# png_of_alphaxTx = pkg_resources.resource_filename('picasso.src.picasso.data', 'gui_data/of_alphaxTx.png')
# png_of_condN = pkg_resources.resource_filename('picasso.src.picasso.data', 'gui_data/of_condN.png')
# png_of_ePe = pkg_resources.resource_filename('picasso.src.picasso.data', 'gui_data/of_ePe.png')
# png_of_ePealphaxTx = pkg_resources.resource_filename('picasso.src.picasso.data', 'gui_data/of_ePealphaxTx.png')
# png_of_sigmaxTsigmax = pkg_resources.resource_filename('picasso.src.picasso.data', 'gui_data/of_sigmaxTsigmax.png')
# png_of_xTx = pkg_resources.resource_filename('picasso.src.picasso.data', 'gui_data/of_xTx.png') 

pyuic5 picasso_design.ui -o picasso_design.py

sed -i "/class Ui_MainWindow/i# inserted via script (begin)" picasso_design.py
sed -i "/class Ui_MainWindow/iimport pkg_resources" picasso_design.py
sed -i "/class Ui_MainWindow/ipng_of_alpha = pkg_resources.resource_filename('picasso.src.picasso.data', 'gui_data/of_alpha.png')" picasso_design.py
sed -i "/class Ui_MainWindow/ipng_of_alphaxTx = pkg_resources.resource_filename('picasso.src.picasso.data', 'gui_data/of_alphaxTx.png')" picasso_design.py
sed -i "/class Ui_MainWindow/ipng_of_condN = pkg_resources.resource_filename('picasso.src.picasso.data', 'gui_data/of_condN.png')" picasso_design.py
sed -i "/class Ui_MainWindow/ipng_of_ePe = pkg_resources.resource_filename('picasso.src.picasso.data', 'gui_data/of_ePe.png')" picasso_design.py
sed -i "/class Ui_MainWindow/ipng_of_ePealphaxTx = pkg_resources.resource_filename('picasso.src.picasso.data', 'gui_data/of_ePealphaxTx.png')" picasso_design.py
sed -i "/class Ui_MainWindow/ipng_of_sigmaxTsigmax = pkg_resources.resource_filename('picasso.src.picasso.data', 'gui_data/of_sigmaxTsigmax.png')" picasso_design.py
sed -i "/class Ui_MainWindow/ipng_of_xTx = pkg_resources.resource_filename('picasso.src.picasso.data', 'gui_data/of_xTx.png')" picasso_design.py
sed -i "/class Ui_MainWindow/i# inserted via script (end)" picasso_design.py
sed -i 's/"of_alpha.png"/png_of_alpha/g' picasso_design.py
sed -i 's/"of_alphaxTx.png"/png_of_alphaxTx/g' picasso_design.py
sed -i 's/"of_condN.png"/png_of_condN/g' picasso_design.py
sed -i 's/"of_ePe.png"/png_of_ePe/g' picasso_design.py
sed -i 's/"of_ePealphaxTx.png"/png_of_ePealphaxTx/g' picasso_design.py
sed -i 's/"of_sigmaxTsigmax.png"/png_of_sigmaxTsigmax/g' picasso_design.py
sed -i 's/"of_xTx.png"/png_of_xTx/g' picasso_design.py

mv picasso_design.py ../../gui/picasso_design.py