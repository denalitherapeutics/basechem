import os

from django.conf import settings
from django.contrib.staticfiles.storage import StaticFilesStorage
from django.core.files.storage import FileSystemStorage

from .storage_backends import MediaStorage, StaticStorage


def select_media_storage():
    if settings.AWS_STORAGE_BUCKET_NAME is None:
        return FileSystemStorage()
    else:
        return MediaStorage()


def select_static_storage():
    if settings.AWS_STORAGE_BUCKET_NAME is None:
        return StaticFilesStorage(location=settings.PROJECT_ROOT + "/public/static")
    else:
        return StaticStorage()


def save_media_file(filepath, fp):
    """
    Saves a file to media storage
    :param filepath: a string, the filepath to save to (excluding MEDIA_ROOT)
    :param fp: a file pointer to the opened file
    :return: a string, the url path to file
    """
    media_storage = select_media_storage()
    if type(media_storage) == FileSystemStorage:
        # If using local storage, make sure path exists
        directory, _ = os.path.split(filepath)
        os.makedirs(f"{settings.MEDIA_ROOT}/{directory}", exist_ok=True)

        # FileSystemStorage does not overwrite existing files but appends a random
        # string to the filepath. To match the behavior with MediaStorage, check and
        # delete the file at the provided path if one exists.
        full_path = os.path.join(settings.MEDIA_ROOT, filepath)
        if os.path.exists(full_path):
            os.remove(full_path)

    media_storage.save(filepath, fp)
    return media_storage.url(filepath)


def upload_to_media_storage(local_filepath):
    """
    Uploads a locally stored file to media storage.
    :param local_filepath: a string, the full path to the file to be uploaded
    """
    media_storage = select_media_storage()
    if not isinstance(media_storage, FileSystemStorage):
        with open(local_filepath, "rb") as file:
            relative_path = os.path.relpath(local_filepath, settings.MEDIA_ROOT)
            media_storage.save(relative_path, file)


def media_file_exists(local_filepath):
    """
    Checks if a file exists on media storage. Downloads it to local storage if it
    doesn't already exist locally.
    :param local_filepath: a string, the full path to the file to be uploaded
    :return: a bool, True if the file exists on media storage
    """
    media_storage = select_media_storage()
    relative_path = os.path.relpath(local_filepath, settings.MEDIA_ROOT)
    exists_on_media_storage = media_storage.exists(relative_path)

    # If the file does not exist locally, open it from storage and save it.
    if not os.path.exists(local_filepath) and exists_on_media_storage:
        os.makedirs(os.path.dirname(local_filepath), exist_ok=True)
        with open(local_filepath, "wb") as file:
            file.write(media_storage.open(relative_path).read())

    return exists_on_media_storage
