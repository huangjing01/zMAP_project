# -*- coding:utf-8 -*-


from flask import Blueprint


reverse_zMap_cancer = Blueprint('reverse_zMap_cancer', __name__)
from . import views
