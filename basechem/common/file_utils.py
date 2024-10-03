def get_tmp_filepath(file, tmp_filename):
    """
    Copies the contents of the given file into a tmp directory with the `tmp_filename`
    :param file: a File object
    :param tmp_filename: a name to give the temporary file (including extension)
    :returns: the path of the tmp file
    """
    tmp_filepath = f"/tmp/{tmp_filename}"
    with open(tmp_filepath, "wb+") as destination:
        for chunk in file.chunks():
            destination.write(chunk)

    return tmp_filepath


def get_tmp_file(f, tmp_filename):
    """
    Adds the given stream to a tmp file with the given filename
    :param f: bytes to write to the file
    :param tmp_filename: a name to give the temporary file (including extension)
    :returns: the path of the tmp file
    """
    tmp_filepath = f"/tmp/{tmp_filename}"
    with open(tmp_filepath, "wb+") as destination:
        destination.write(f)

    return tmp_filepath
