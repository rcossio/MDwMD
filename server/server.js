const express = require('express');
const bodyParser = require('body-parser');
const fs = require('fs');
const path = require('path');

const app = express();
const PORT = 3000;

// Middleware
app.use(bodyParser.json());
app.use(express.static('public'));

// Serve HTML Page
app.get('/', (req, res) => {
  res.sendFile(path.join(__dirname, '/public/index.html'));
});

// Handle POST Request
app.post('/update-data', (req, res) => {
  const newData = req.body;
  const filePath = path.join(__dirname, 'loaded_data.json');

  fs.readFile(filePath, (err, data) => {
    if (err) {
      res.status(500).send('Error reading data file');
      return;
    }

    const jsonData = JSON.parse(data);
    jsonData.push(newData);

    fs.writeFile(filePath, JSON.stringify(jsonData, null, 2), (err) => {
      if (err) {
        res.status(500).send('Error updating data file');
        return;
      }

      res.send({status:'success',payload:'Data updated successfully'});
    });
  });
});

app.listen(PORT, () => {
  console.log(`Server is running on http://localhost:${PORT}`);
});
