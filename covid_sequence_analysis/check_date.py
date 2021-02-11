import datetime

def check_date_format(input_date):

    """ Check the format of a date to determine if it matches YYYY-MM-DD. Returns a tuple containing a bool
    (conditional on whether or not the date is in a valid format) and, if True, the formatted date as a str, 
    or if False a ValueError. """

    is_valid_date = True
    try:
        year, month, day = input_date.split('-')   
        datetime.datetime(int(year),int(month),int(day))
    except ValueError:
        is_valid_date = False

    if is_valid_date:
        return is_valid_date, str(datetime.datetime(int(year),int(month),int(day)).date())
    else :
        return is_valid_date, ValueError('This is not a valid date format. Use the format YYYY-MM-DD.')


