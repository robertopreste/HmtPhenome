#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
# from flask_wtf import Form, FlaskForm
# from wtforms import StringField, SelectField, SelectMultipleField, RadioField, BooleanField, PasswordField, widgets
# from wtforms.validators import InputRequired, Length, NumberRange, Optional, AnyOf, DataRequired, Email, EqualTo, ValidationError
# from .models import User


# class MultiCheckboxField(SelectMultipleField):
#     widget = widgets.ListWidget(prefix_label=False)
#     option_widget = widgets.CheckboxInput()


# class QueryForm(Form):
#     first_name = StringField("first_name", validators=[DataRequired()])
#     last_name = StringField("last_name", validators=[DataRequired()])
#     country = StringField("country", validators=[DataRequired()])
#
#
# class RegistrationForm(FlaskForm):
#     username = StringField("username", validators=[DataRequired()])
#     email = StringField("email", validators=[DataRequired(), Email()])
#     password = PasswordField("password", validators=[DataRequired()])
#     password2 = PasswordField("password2", validators=[DataRequired(), EqualTo("password")])
#     first_name = StringField("first_name", validators=[DataRequired()])
#     last_name = StringField("last_name", validators=[DataRequired()])
#
#     def validate_username(self, username):
#         user = User.query.filter(User.username == username.data).first()
#         if user is not None:
#             raise ValidationError("Please choose a different username.")
#
#     def validate_email(self, email):
#         user = User.query.filter(User.email == email.data).first()
#         if user is not None:
#             raise ValidationError("Please choose a different email address.")
#
#
# class LoginForm(FlaskForm):
#     username = StringField("username", validators=[DataRequired()])
#     password = PasswordField("password", validators=[DataRequired()])


