import factory
from django.core.files.uploadedfile import SimpleUploadedFile
from django_q.models import OrmQ
from rdkit import Chem

from basechem.main.models.collection_models import Collection
from basechem.main.models.compound_models import Compound, CompoundOccurrence, Series
from basechem.main.models.project_models import Project
from basechem.users.tests.factories import BasechemUserFactory


class ProjectFactory(factory.DjangoModelFactory):
    code = factory.Faker("lexify", text="????")

    class Meta:
        model = Project


class CollectionFactory(factory.DjangoModelFactory):
    owner = factory.SubFactory(BasechemUserFactory)
    project = factory.SubFactory(ProjectFactory)

    class Meta:
        model = Collection


class CompoundFactory(factory.DjangoModelFactory):
    class Meta:
        model = Compound


class CompoundOccurrenceFactory(factory.DjangoModelFactory):
    compound = factory.SubFactory(CompoundFactory)
    owner = factory.SubFactory(BasechemUserFactory)
    parent_co = None

    @factory.post_generation
    def parent_co_update(self, create, extracted, **kwargs):
        """
        If `parent_co` was specified,  set the `compound` and `owner` to match the `parent_co`
        """
        if self.parent_co:
            self.owner = self.parent_co.owner
            self.compound = self.parent_co.compound
            self.save()

    class Meta:
        model = CompoundOccurrence


class SeriesFactory(factory.DjangoModelFactory):
    project = factory.SubFactory(ProjectFactory)

    class Meta:
        model = Series

    @classmethod
    def _create(cls, model_class, *args, **kwargs):
        """
        Add a bound state file from smiles if one doesn't exist already
        """
        obj = model_class(*args, **kwargs)

        if not obj.bound_state_file:
            mol = Chem.MolFromSmiles(obj.smiles)
            tmp_path = "/tmp/series_test.sdf"
            writer = Chem.SDWriter(tmp_path)
            writer.write(mol)
            writer.close()

            obj.bound_state_file = SimpleUploadedFile(
                f"/tmp/{obj.dn_id}_series.sdf", open(tmp_path, "rb").read()
            )

        obj.save()
        return obj


class OrmQFactory(factory.DjangoModelFactory):
    class Meta:
        model = OrmQ
