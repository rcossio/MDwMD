const express = require('express');
const bodyParser = require('body-parser');
const mongoose = require('mongoose');
const xml2js = require('xml2js');
const Bottleneck = require('bottleneck');
require('dotenv').config();

const app = express();
const PORT = 3000;

// Initialize a new bottleneck limiter
const limiter = new Bottleneck({
  maxConcurrent: 1, // Maximum number of jobs running at the same time
  minTime: 500 // Minimum time (in milliseconds) between job executions
});

// MongoDB Schema and Model setup
const diffusionSchema = new mongoose.Schema({
  diffusionCoefficient: Number,
  diffusionError: Number,
  diffusionUnit: String,
  manualProcessing: String,
  otherManualProcessing: String,  
  method: String,
  otherMethod: String,
  proteinName: String,
  species: String,
  accessionNumber: String,
  referenceId: String,
  referenceIdType: String,
  referenceType: String,
  comment: String, 
  timestamp: String,
  userName: String,
  active: Boolean,
  supersededBy: { type: mongoose.Schema.Types.ObjectId }
});

// Create Mongoose Model
const DiffusionData = mongoose.model('data_from_experiments', diffusionSchema);

app.use(bodyParser.json());
app.use(express.static('public'));

// Serve HTML Page
app.get('/', (req, res) => {
  res.sendFile(path.join(__dirname, '/public/index.html'));
});

// Handle POST Request
app.post('/update-data', async (req, res) => {
  const {data, replacedEntry} = req.body;

  try {
    const newEntryId = await DiffusionData.create({...data, active: true}); 
    if(replacedEntry){
      await DiffusionData.updateOne({_id: replacedEntry}, { $set: { supersededBy: newEntryId, active: false } });
    }
    res.send({ status: 'success', payload: 'Data saved to MongoDB' });
  } catch (err) {
      console.error('Error saving data:', err);
      res.status(500).send('Error saving data');
  }
});

// Handle GET Request for Search
app.get('/search', (req, res) => {
  const { query,species } = req.query;

  let queryConditions = [
    { active: true },
    {
      $or: [
        { proteinName: new RegExp(query, 'i') },
        { accessionNumber: new RegExp(query, 'i') }
      ]
    }
  ];
  
  // Add species filter only if speciesList is not empty and not null
  const speciesList = species.split(',').map(item => item.trim()).filter(item => item !== '')
  
  if (speciesList && speciesList.length > 0) {
    queryConditions.push({ species: { $in: speciesList } });
  }


  DiffusionData.find({
    $and: queryConditions
  })
  .sort({ accessionNumber: 1 })
  .then(data => {
    // Extract species from the results and create a unique list
    const species = new Set(data.map(item => item.species));
    const speciesList = Array.from(species).map(Number).sort((a, b) => a - b);    

    // Send the response with both results and species list
    res.json({
      status: 'success',
      payload: {
        result: data,
        species: speciesList
      }
    });
  })
  .catch(err => {
    console.error('Error during search:', err);
    res.status(500).send('Error during search');
  });
});

// Handle GET Request for fetching publication year
app.get('/publication-year', (req, res) => {
  const { refId, refType } = req.query;

  const fetchPublicationYear = async () => {
    let url, response, text;
    if (refType === 'pmid') {
      url = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=${refId}&retmode=xml`;
    } else if (refType === 'doi') {
      url = `https://api.crossref.org/works/${refId}`;
    }
    
    try {
      const fetch = (await import('node-fetch')).default;
      response = await fetch(url);
      text = await response.text();
    } catch (error) {
      console.error('Fetch error:', error);
      return Promise.reject('Failed to fetch data');
    }

    if (refType === 'pmid') {
      return new Promise((resolve, reject) => {
        xml2js.parseString(text, (err, result) => {
          if (err) {
            reject('Error parsing XML');
          } else {
            try {
              const year = result.PubmedArticleSet.PubmedArticle[0].MedlineCitation[0].Article[0].Journal[0].JournalIssue[0].PubDate[0].Year[0];
              resolve(year);
            } catch (parseError) {
              reject('Error extracting publication year');
            }
          }
        });
      });
    } else if (refType === 'doi') {
      try {
        const json = JSON.parse(text);
        const year = json.message['published-print']?.['date-parts'][0][0] || json.message['published-online']?.['date-parts'][0][0];
        return Promise.resolve(year);
      } catch (parseError) {
        console.error('Error parsing JSON:', parseError);
        return Promise.reject('Error extracting publication year from DOI');
      }
    }
  };

  limiter.schedule(fetchPublicationYear)
    .then(year => res.json({ year }))
    .catch(error => {
      console.error(error);
      res.status(500).send(error);
    });
});

mongoose.connect(process.env.MONGODB_URI)
.then(() => console.log('Connected to MongoDB'))
.catch(err => console.error('Error connecting to MongoDB:', err));

app.listen(PORT, () => {
  console.log(`Server running on http://localhost:${PORT}`);
});
