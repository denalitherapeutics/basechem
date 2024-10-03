# pragma: no cover

from django.contrib.postgres.operations import CreateExtension
from django.db import models


class RDKitExtension(CreateExtension):
    def __init__(self):
        self.name = "rdkit"

    def database_forwards(self, app_label, schema_editor, from_state, to_state):
        super(RDKitExtension, self).database_forwards(
            app_label, schema_editor, from_state, to_state
        )
        # Ensure the RDKit extension is fully loaded by running command
        with schema_editor.connection.cursor() as c:
            c.execute("SELECT mol_to_smiles('C'::mol)")


class MolField(models.Field):
    description = "A rdkit mol object"

    def __init__(self, *args, **kwargs):
        kwargs["help_text"] = "Field to manage RDKit Mol Objects"
        super().__init__(*args, **kwargs)

    def db_type(self, connection):
        return "mol"
