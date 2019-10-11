#!/usr/bin/env python3
"""
Contains functions for grabbing data from the Internet.

TUNL:
- we download the necessary pdf from TUNL
- then look through the PDF and pick out important values using tabula

But TUNL wasn't a great source since their data isn't very parse-able

Instead, use Brookhaven's site: nndc.bnl.gov 


"""
from os import mkdir
from os.path import join, realpath, split, exists
from .ncsd_output_reader import element_name

this_dir = split(realpath(__file__))[0]
    
def parse_tunl_pdf(pdf_url, pdf_save_path):
    """returns Ex column of pdf from TUNL, without uncertainties"""
    # I've done something gross here by putting some imports inside functions
    # rather than up here... The reason for that is that some packages,
    # especially tabula, are not always going to be installed.
    # This way, if the user wants to install them they can,
    # but if not, the code in this package still works.
    import urllib.request
    import requests
    from lxml import html
    import tabula
    print("\nThe feature to grab data from TUNL is still very experimental,\n"+\
        "please ensure you end up with the right values!\n"+\
        "Also be sure to manually adjust J and T values.\n\n")
    print("downloading PDF")
    # saves url file to a local destination
    urllib.request.urlretrieve(pdf_url, pdf_save_path)
    print("pdf saved as "+pdf_save_path)
    print("reading into dataframe")
    # read in the PDF file that contains data
    df = tabula.read_pdf(pdf_save_path, pages="all")
    # I tried a couple things but wasn't sure how to get the J or T values...
    # get the relevent column
    try:
        Ex_strings = df["Ex (MeV± keV)"]
    except KeyError:
        Ex_strings = df["E x"]
    Ex_strings.dropna(inplace = True)  # drop all NaN entries
    
    Exs = []  # to store the Ex values
    for Ex_string in list(Ex_strings):
        # there are so many things that could go wrong with this...
        # but let's try anyway!
        
        # first off, sometimes we just get random integers.
        # let's ignore those by ensuring all lines have a decimal point
        # except for sometimes the first value is zero and I assume we want that
        if "." not in Ex_string and Ex_string != "0":
            continue
        
        try:
            Ex = float(Ex_string)
            Exs.append(Ex)
            continue
        except ValueError:
            pass
        # if we're still here, 
        # we could not convert string to float for some reason...

        # there are sometimes artifacts in the data from super/subscripts
        # thankfully those always seem to come with + or - signs, e.g. 1+
        # and they're just 1 character as far as I've seen
        for abomination in ["-", "+"]:
            while abomination in Ex_string:
                index = Ex_string.index(abomination)
                # I assume none of these + or - occur right at the start
                before = Ex_string[:index-1]  # also ignore the char before
                after = Ex_string[index+1:]
                Ex_string = before + after
        try:
            Ex = float(Ex_string)
            Exs.append(Ex)
            continue
        except ValueError:
            pass    

        # remove uncertainties from the string
        if "±" in Ex_string:
            index = Ex_string.index("±")
            before = Ex_string[:index]
            after = Ex_string[index+1:]
            Ex_string = before  # just take the stuff before the ±
        try:
            Ex = float(Ex_string)
            Exs.append(Ex)
            continue
        except ValueError:
            pass

        # if it STILL doesn't work, just remove all nonsensical characters
        # and hope that works
        new_string = ""
        for char in Ex_string:
            if (not char.isdigit()) and (char != "."):
                pass
            else:
                new_string += char
        if new_string != "":
            Ex = float(new_string)
            Exs.append(Ex)
    return Exs


def get_tunl_data(calc_data):
    """Function to get data from TUNL. Currently only gets energies,
    but if you have any ideas about how to get J, T, parity let me know!"""
    import urllib.request
    import requests
    from lxml import html
    import tabula
    
    # a place to save all the pdfs we download, we'll overwrite for new files
    save_dir = realpath(join(this_dir, '..'))
    pdf_save_path = join(save_dir, 'TUNL.pdf')
    
    Z = calc_data["Z"]
    N = calc_data["N"]
    n_states = calc_data["n_states"]

    A = Z + N
    element = element_name(Z)
    print("getting TUNL data for "+element+str(A))
    # get page for elements with our A value
    page = requests.get(
        "http://www.tunl.duke.edu/nucldata/chain/{A}.shtml".format(A=A))
    tree = html.fromstring(page.content)

    # find all elements with the right element name, e.g. "Li" or "B"
    elements = tree.xpath("//*[text()='{element}']".format(element=element))

    # find all possible pdf files associated with this element
    possible_files = []
    for element in elements:
        link = element.attrib.get('href')
        if ".pdf" in link:
            possible_files.append(link)
    
    # check that we only got one
    if len(possible_files) != 1:
        print("Possible PDF files found on this page:")
        print(possible_files)
        raise ValueError(
            "Wrong number of files found on TUNL: "+str(len(possible_files)))
    
    # then take that file, open it, and look inside for values
    pdf_url = possible_files[0]
    Ex_values = parse_tunl_pdf(pdf_url, pdf_save_path)  # list of floats
    
    # format the data for use later
    data = {
        "expt_spectrum": {
            "Expt": {
                n+1: [1,1,1, Ex_values[n]]
                for n in range(n_states) # assuming we want the first few
            }       
        }       
    }
    return data


def get_bnl_data(calc_data):
    """Function to get data from BNL"""
    from selenium import webdriver
    from selenium.webdriver.common.keys import Keys
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait
    from selenium.webdriver.support import expected_conditions as EC

    Z = calc_data["Z"]
    N = calc_data["N"]

    A = Z + N
    element = element_name(Z)
    print("getting BNL data for "+element+str(A))
    path_to_chromedriver = join(this_dir, "chromedriver")
    driver = webdriver.Chrome(path_to_chromedriver)
    # open main search page
    driver.get("https://www.nndc.bnl.gov/ensdf/")
    # get page for elements with our A value
    search_box = driver.find_element_by_name("nuc")
    search_box.clear()
    search_box.send_keys(str(A))
    search_box.send_keys(Keys.RETURN)
    
    # find the right box to click
    rows = (driver.find_elements_by_class_name("normrow-ensdf") +
            driver.find_elements_by_class_name("altrow-ensdf"))
    for row in rows:
        html_text = row.get_attribute('innerHTML')
        words = html_text.split()
        if element in words and "ADOPTED LEVELS" in html_text:
            checkbox = row.find_element_by_name("datasetcheck")
            checkbox.click()
            break
    # click button to see ENSDF text format
    xpath = "//input[@name='chooseit' and @value='ENSDF text format']"
    text_button = driver.find_elements_by_xpath(xpath)[0]
    text_button.click()

    # switch tabs (to the 1th tab, i.e. the new one that just opened)
    driver.switch_to.window(driver.window_handles[1])

    # now we have the page with text, get text
    text = driver.find_element_by_tag_name("body").text

    # split into lines, grab data
    prefix = str(A) + element
    lines = text.split("\n")
    online_data = []
    for line in lines:
        line = line[len(prefix)+2:]  # ignore prefix
        line_info = line[:3]
        
        # first line_info[0] says what kind of line it will be
        if line_info[0] == " ":
            line_type = 1
        elif line_info[0] == "2" and "ISPIN" in line:
            line_type = 2
        else:
            continue
        # line_info[0] designates comments
        if line_info[1] == "c":
            continue
        # line_info[2] must be "L" for it to matter
        if line_info[2] != "L":
            continue
        
        # if line contains any greater than or less than symbols remove them
        line = line.replace("GE", "").replace("LE", "")

        line = line[len(line_info):]  # take only the rest of the line
        words = line.split()
        if len(words) <= 1:
            continue

        # we have two options, either isospin line or energy line
        if line_type == 1:
            energy = float(words[0])
            if "+" not in words[1] and "-" not in words[1]:
                # we have two options: either there's an extra word here
                # or we don't have parity
                word2 = words[2]
                if "+" in word2 or "-" in word2:
                    word = word2
                else:
                    word = words[1]
            else:
                word = words[1]
            word = word.replace("(", "").replace(")", "")
            if "," in word:
                word = word[:word.index(",")]
            # try to get parity
            parity = word[-1]
            if parity in ["+", "-"]:
                if parity == "+":
                    parity = 1
                else:
                    parity = -1
                word = word[:-1] # chop off parity before trying to get J
            else:
                parity = "?"
            # now get J, dealing with fractions if necessarity
            if "/" in word:
                num, den = word.split("/")
                word = str(float(num) / float(den))
            j = float(word)
        else:
            word = words[0][6:]
            word = word.replace("(", "").replace(")", "")
            for trash in ["$", "+", "-"]:
                if trash in word:  # remove it
                    word = word[:word.index(trash)]
            if "/" in word:  # case for if there's a number of the form 1/2
                num, den = word.split("/")
                word = str(float(num) / float(den))
            isospin = float(word)
            online_data.append([j, isospin, parity, energy])

    calc_spectrum = calc_data["calculated_spectrum"]

    # output data format
    online_dict = {"expt_spectrum": { "Expt": { } } }
    # in Expt the format will be:
    # state_num: [2J, 2T, parity, energy]

    # produce experimental spectrum
    nmax_max = max(calc_spectrum.keys())
    for state in sorted(calc_spectrum[nmax_max].keys()):
        J2, repetition, parity, energy = calc_spectrum[nmax_max][state]
        for o_d in online_data:
            # o_ for "online"
            o_J, o_T, o_parity, o_energy = o_d
            o_J2 = int(2 * o_J)
            o_T2 = int(2 * o_T)
            # take online energy if online data matches all other parameters
            if J2 == o_J2 and parity == o_parity:
                energy = o_energy
                break
        online_dict["expt_spectrum"]["Expt"][state] = [
            J2, repetition, parity, energy]

    return online_dict


def get_online_data_wrapper(calc_data):
    """This tries to get online data,
    but doesn't break everything if it fails.
    Instead, filler data is returned.
    """
    src = "bnl"
    name_map = {"bnl": "BROOKHAVEN", "tunl": "TUNL"}
    func_map = {"bnl": get_bnl_data, "tunl": get_tunl_data}
    try:
        print(1/"a")
        data = func_map[src](calc_data)
        return data
    except Exception as e:
        print("Warning raised when trying to get data from "+name_map[src]+":")
        print(e)
        print("if your exception was 'no module named <x>', "
              "try runnning 'pip install --user <x>'")
        print("Returning filler data instead")
        # produce experimental spectrum
        calc_spectrum = calc_data["calculated_spectrum"]
        nmax_max = max(calc_spectrum.keys())
        filler_data = { "expt_spectrum": { "Expt": { } } }
        for state in sorted(calc_spectrum[nmax_max].keys()):
            filler_data["expt_spectrum"]["Expt"][state] = \
                calc_spectrum[nmax_max][state]
        return filler_data
